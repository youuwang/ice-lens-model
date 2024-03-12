#include "Ice_model.hpp"
#include "FDM_ff.hpp"
#include "FDM_icelens.hpp"
#include "input_file.hpp"

using namespace std;
class Ice_model;
void Ice_model::read_input(std::string filename) {
    // define the section you are looking for: name it "section"
    const std::vector<std::string> sections = { "material", "nodes", "domain", "initial_condition", "BC" };

    // read lines (all sections) from ifile and store in container name "lines"
    // std::vector<std::string> lines = read_input_section(filename, section);
    std::vector<std::string> lines;
    for (auto section : sections) {
        auto lns = read_input_section(filename, section);
        lines.insert(lines.end(), lns.begin(), lns.end());
        }

        // loop over content of lines and find corresponding variables
        for (std::string line: lines) {
        // std::cout << line << std::endl;
        std::stringstream sstr(line);
        std::string word;

        sstr >> word;
        if (word == "$mat_prop") {
            sstr >> this->rho_s;
            sstr >> this->cp_s;
            sstr >> this->phi;
            sstr >> this->c;
            sstr >> this->k0;
            sstr >> this->Rp;
        } 
        else if (word == "$nodes") {
            sstr >> this->N;
        } 
        else if (word == "$length") {
            sstr >> this->d;
        } 
        else if (word == "$temperature") {
            sstr >> this->T01;
            sstr >> this->T02;
        }
        else if (word == "$bc") {
            sstr >> this->T1;
            sstr >> this->T2;
            sstr >> this->v_ch;
            sstr >> this->cyc_tot;
        }
    }
}


void Ice_model::solve(ofstream& output_T, ofstream& output_ice_lens, ofstream& label_time, ofstream& crack_state_time) {
    ofstream time_N_lens("t_N.dat");
    ofstream ub_lb("ub_lb_t");
    ofstream time_xf2("t_xf2.dat");
    ofstream time_xf3("t_xf3.dat");
    ofstream Tu_Tl("Tu_Tl.dat");
    ofstream time1("t1.dat");
    ofstream time_dry("t_dry.dat");
    ofstream active_ice_lens("active_ice_lens.dat");

    while(this->N_cyc < this->cyc_tot) {
        switch (this->flag) {
            case 0:
                if (this->N_cyc == 1) {
                    this->solve3_active = false;
                    this->solve1(output_T);
                    std::cout << "Phase 1 of cycle " << this->N_cyc << " finished ..." << std::endl;
                    if (this->time > this->v_ch / 2) {
                        std::cout << "No frozen fringe formed " << std::endl;
                        this->flag = 1;
                        goto case1;
                    }

                    this->solve2(output_T);
                    if (this->time > this->v_ch / 2) {
                        std::cout << "Frozen fringe involved but no ice lens formed in cycle" << this->N_cyc << std::endl;
                        this->flag = 1;
                        goto case1;
                    }

                    this->solve3(output_T, label_time, crack_state_time);
                    this->flag = 1;
                }
                else {
                    this->solve3(output_T, label_time, crack_state_time);
                    this->flag = 1;
                    std::cout << "Cooling ends" << std::endl;
                }
            
            case 1:
            case1:
                if (this->solve3_active == true) {
                    this->solve3(output_T, label_time, crack_state_time);
                }
                else if (this->solve2_active == true) {
                    this->solve2(output_T);
                }
                else {
                    this->solve1(output_T);
                }
                this->flag = 0;
        }

        std::cout << "Cycle " << this->N_cyc << " finished and current ice lenses are:" << std::endl;
        for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
            std::cout << this->x_ub[i] << "\t" << this->x_lb[i] << std::endl;
        }
        this->N_cyc += 1;
    }

    //! output results
    for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
        output_ice_lens << this->x_ub[i] << "\t" << this->x_lb[i] << endl;
    }

    if (this->time_t2.size() > 0) {
        for (int m = 0; m < static_cast<int>(this->xf_t2.size()); m++) {
            time_xf2 << this->time_t2[m] << "\t" << this->xf_t2[m] << endl;
        }
    }

    for (int n = 0; n < static_cast<int>(this->xf_t3.size()); n++) {
        time_xf3 << this->time_t3[n] << "\t" << this->xf_t3[n] << endl;
    }

    for (int j = 0; j < static_cast<int>(this->time_t3.size()); j++) {
        time_N_lens << this->time_t3[j] << "\t" << this->N_lens[j] << endl;
    }

    for (int k = 0; k < static_cast<int>(this->x_ub_t.size()); k++) {
        ub_lb << this->x_ub_t[k] << "\t" << this->x_lb_t[k] << endl;
    }

    for (int l = 0; l < static_cast<int>(this->T_ub_t.size()); l++) {
        Tu_Tl << this->T_ub_t[l] << "\t" << this->Tl_t[l] << endl;
    }

    for (int p = 0; p < static_cast<int>(this->time_t1.size()); p++) {
        time1 << this->time_t1[p] << endl;
    }

    for (int q = 0; q < static_cast<int>(this->x_dry_t.size()); q++) {
        time_dry << this->time_t3[q] << "\t" << this->x_dry_t[q] << endl;
    }

    for (int r = 0; r < static_cast<int>(this->time_t3.size()); r++) {
        active_ice_lens << this->time_t3[r] << "\t" << this->ice_active_t[r] << endl;
    }

    time_N_lens.close();
    ub_lb.close();
    time_xf2.close();
    time_xf3.close();
    Tu_Tl.close();
    time1.close();
    time_dry.close();
    active_ice_lens.close();
}


void Ice_model::solve1(ofstream& output_T) {

    this->solver_ptr = std::make_shared<FDM_solver>(this->d, this->N);

    // Solve the heat equation using the finite difference method
    while (true) {
        this->time_temp += this->dt;
        this->time += this->dt;

        if (this->T[0] < (this->Tf - 273.15) || (this->T[0] > this->T2)) {
            this->solve2_active = true;
            return;
        }

        if (this->time > this->v_ch / 2) {
            return;
        }

        //output results
        if (fabs(this->time_output - 0.1) < 1e-8) {
            for (int i = 0; i < this->N; i++) {
                if (i < this->N - 1) {
                    output_T << this->T[i] << "\t";
                }
                else {
                    output_T << this->T[i] << endl;
                }
            }

            this->time_t1.push_back(this->time);
            this->time_output = 0;
        }
        this->time_output += this->dt;

        // // bi-linear temperature
        // if (this->flag == 0) {
        //     this->T[0] -= this->v_ch * this->dt;
        // }
        // else {
        //     this->T[0] += this->v_ch * this->dt;
        // }

        // Cosine temperature
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_ch);

        this->solver_ptr->solve(this);
    }
}

void Ice_model::solve2(ofstream& output_T) {
    double integral; // integral inside the expression of Fp
    int N_temp = 100;
    std::vector<double> Fp(N_temp, 0); // net interparticle force at each node
    std::vector<double> x_temp(N_temp, 0); // temporary x for calculation of Fp
    int index_x_temp; // index of element where x_temp is located
    double Tx; // temperature at x_temp
    std::vector<double> T_temp;

    this->solver_ptr = std::make_shared<FDM_ff>(this->d, this->N);

    ofstream F_p("F_p.dat");

    // Solve the heat equation using the finite difference method
    while (true) {
        this->time += this->dt;
        this->time_temp += this->dt;

        // std::vector<double> x0_temp;
        // std::vector<int> idx_0_temp;

        // Criterion for ending the second phase
        if (this->time_temp > this->v_ch / 2) {
            this->time_temp = 0;
            return;
        }

        T_temp = this->T;

        //locate Tm
        for (int i = 0; i < this->N - 1; i++) {
            if ((T_temp[i] * T_temp[i + 1] < 0) || (T_temp[i] == 0)) {
                this->x0_active = this->x[i] + this->dx * (- T_temp[i]) / (T_temp[i + 1] - T_temp[i]);
            }
        }

        // locate Tf
        for (int i = 0; i < this->N - 1; i++) {
            if (((this->Tf - 273.15 - this->T[i]) * (this->Tf - 273.15 - this->T[i + 1]) < 0) || (this->Tf - 273.15 - this->T[i] == 0)) {
                this->xf = this->x[i] + this->dx * ((this->Tf - 273.15) - this->T[i]) / (this->T[i + 1] - this->T[i]);
            }
        }

        this->idx = this->find_index(this->xf);

        // get net interparticle force profile above frozen fringe boundary
        for (int i = 0; i < N_temp; i++) {
            x_temp[i] = this->xf / (N_temp - 1) * i;
        }
        for (int j = 0; j < N_temp; j++) {
            index_x_temp = this->find_index(x_temp[j]);
            integral = this->FT_int(x_temp[j]);
            Tx = this->T[index_x_temp] + (x_temp[j] - this->x[index_x_temp]) / this->dx * (this->T[index_x_temp + 1] - this->T[index_x_temp]) + 273.15;
            Fp[j] = this->c * (1 - this->phi) - this->rho_i * this->L / this->Tm * (integral + (this->Tm - Tx) * (1 - this->phi * this->ice_sat_fun(Tx)));
        }

        // Output the result
        if (fabs(this->time_output - 0.1) < 1e-8) {
            for (int i = 0; i < this->N; i++) {
                if (i < this->N - 1) {
                    output_T << this->T[i] << "\t";
                }
                else {
                    output_T << this->T[i] << endl;
                }
            }
            this->time_t2.push_back(this->time);
            this->xf_t2.push_back(this->xf);
            this->time_output = 0;
            for (int j = 0; j < N_temp; j++) {
                if (j < N_temp - 1) {
                    F_p << Fp[j] << "\t";
                }
                else {
                    F_p << Fp[j] << endl;
                }
            }
        }

        this->time_output += this->dt;

        // // bi-linear T
        // if (this->flag == 0) {
        //     this->T[0] -= this->v_ch * this->dt;
        // }
        // else {
        //     this->T[0] += this->v_ch * this->dt;
        // }

        // Cosine T
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_ch);

        // call the solver
        this->solver_ptr->solve(this);

        for (int j = 0; j < N_temp - 1; j++) {
            if ((Fp[j] * Fp[j + 1] < 0) || (Fp[j] == 0)) {
                this->lens_ini = x_temp[j] + (-Fp[j]) / (Fp[j + 1] - Fp[j]) * this->xf / (N_temp - 1);
                std::cout << "Ice lens formed at " << this->lens_ini << std::endl;
                this->solve3_active = true;
                this->x_lb.push_back(this->lens_ini);
                this->x_ub.push_back(this->lens_ini);
                this->idx_lb.push_back(this->find_index(this->lens_ini));
                this->idx_ub.push_back(this->find_index(this->lens_ini));
                this->crack_state.push_back(0);

                //! calculate the temperature of the first ice lens
                int ice_active_ini = this->find_index(this->lens_ini);
                this->Tl.push_back(this->T[ice_active_ini] + (lens_ini - this->x[ice_active_ini]) / this->dx * (this->T[ice_active_ini + 1] - this->T[ice_active_ini]) + 273.15);
                this->T_ub.push_back(this->T[ice_active_ini] + (lens_ini - this->x[ice_active_ini]) / this->dx * (this->T[ice_active_ini + 1] - this->T[ice_active_ini]) + 273.15);                    

                // locate Tm
                for (int i = 0; i < this->N - 1; i++) {
                    if ((T_temp[i] * T_temp[i + 1] < 0) || (T_temp[i] == 0)) {
                        this->x0.push_back(this->x[i] + this->dx * (- T_temp[i]) / (T_temp[i + 1] - T_temp[i]));
                    }  
                }

                return;
            }
        }     
    }
    F_p.close();
}

void Ice_model::solve3(ofstream& output_T, ofstream& label_time, ofstream& crack_state_time) {
    // std::vector<double> x0_prev; // previous position of 0 degree Celsius
    double integral_FT; // integral inside the expression of FT
    double integral_F_mu; // integral inside the expression of F_mu
    int N_temp = 50; // size of x_temp
    std::vector<double> Fp; // net interparticle force at each node of interest
    double x_temp; // temporary x for calculation of Fp
    int index_x_temp; // index of element where x_temp is located
    double Tx; // temperature at x_temp
    this->solver_ptr = std::make_shared<FDM_icelens>(this->d, this->N);
    std::vector<double> T_temp;
    bool new_ice_lens; // indicator for new ice lens formation
    std::vector<double> x0_temp;
    std::vector<double> mw; // water amount stored
    std::vector<double> x_lb_temp; // temporary x_lb
    std::vector<double> x_ub_temp; // temporary x_ub
    std::vector<double> xf_all; // all xf
    std::vector<double> x_temp1; // vector of x_temp

    while (true) {
        this->time += this->dt;
        this->time_temp += this->dt;

        if (this->time_temp > this->v_ch / 2) {
            this->time_temp = 0;
            return;
        }

        auto xl_max = std::max_element(this->x_lb.begin(), this->x_lb.end());
        if (*xl_max > this->d) {
            std::cout << "Over" << std::endl;
            return;
        }

        for (int i = 0; i < this->N - 1; i++) {
            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                if ((this->x[i] - this->x_ub[j]) * (this->x[i] - this->x_lb[j]) < 0) {
                    this->phi_n[i] = 1;
                    break;
                }
            }
            if (this->phi_n[i] == this->phi) {
                if ((this->x[i] > this->x_dry) || (this->x[i] == this->x_dry)) {
                    this->label[i] = 0;
                }
                else {
                    this->label[i] = 1;
                }
            }
            else if (this->T[i] > 0) {
                this->label[i] = 3;
            }
            else {
                this->label[i] = 2;
            }
        }
        this->phi_n.assign(this->N, this->phi);

        // locate Tm
        for (int i = 0; i < this->N - 1; i++) {
            if ((this->T[i] * this->T[i + 1] < 0) || (this->T[i] == 0)) {
                x0_temp.push_back(this->x[i] + this->dx * (- this->T[i]) / (this->T[i + 1] - this->T[i]));
            }
        }
        this->x0 = x0_temp;
        x0_temp.clear();

        for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
            if ((this->Tl[i] - this->Tm) * (this->T_ub[i] - this->Tm) < 0) {
                this->crack_state[i] = 1;
            }
            else if (this->T_ub[i] > this->Tm) {
                this->crack_state[i] = 2;
            }
            else {
                this->crack_state[i] = 0;
            }
        }

        this->ice_active = this->ice_lens_active();

        // locate Tf
        for (int i = 0; i < this->N - 1; i++) {
            if (((this->Tf - 273.15 - this->T[i]) * (this->Tf - 273.15 - this->T[i + 1]) < 0) || ((this->Tf - 273.15) == this->T[i])) {
                this->xf = this->x[i] + this->dx * ((this->Tf - 273.15) - this->T[i]) / (this->T[i + 1] - this->T[i]);
                this->idx = i;
                xf_all.push_back(this->xf);
            }
        }

        // determine the active x0
        if (this->ice_active != 1e6) {
            for (int j = 0; j < static_cast<int>(this->x0.size()); j++) {
                if (this->x0[j] > this->x_lb[this->ice_active]) {
                    this->x0_active = this->x0[j];
                    break;
                }
            }
        }

        T_temp = this->T;

        if (fabs(this->time_output - 0.1) < 1e-8) {
            for (int i = 0; i < this->N; i++) {
                if (i < this->N - 1) {
                    output_T << T_temp[i] << "\t";
                }
                else {
                    output_T << T_temp[i] << endl;
                }
            }

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                if (j < static_cast<int>(this->x_lb.size()) - 1) {
                    crack_state_time << this->crack_state[j] << "\t";
                }
                else {
                    crack_state_time << this->crack_state[j] << endl;
                }
            }

            for (int j = 0; j < static_cast<int>(this->x_ub.size()); j++) {
                this->x_ub_t.push_back(this->x_ub[j]);
                this->x_lb_t.push_back(this->x_lb[j]);
                this->T_ub_t.push_back(this->T_ub[j]);
                this->Tl_t.push_back(this->Tl[j]);
            }

            this->time_t3.push_back(this->time);
            this->xf_t3.push_back(this->xf);
            this->x_dry_t.push_back(this->x_dry);
            this->ice_active_t.push_back(this->ice_active);

            this->N_lens.push_back(static_cast<int>(this->x_ub.size()));

            if (this->water_supply == 0) {
                std::cout << "No water supply" << std::endl;
            }

            for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
                if (this->ice_active == 1e6) {
                    std::cout << this->T_ub[i] << "\t" << this->Tl[i] << "\t" << "No active ice lenses" << std::endl;
                }
                else {
                    std::cout << this->T_ub[i] <<"\t"<< this->Tl[i] << "\t" << this->vl << std::endl;
                }
            }

            for (int i = 0; i < this->N; i++) {
                if (i < this->N - 1) {
                    label_time << this->label[i] << "\t";
                }
                else {
                    label_time << this->label[i] << endl;
                }
            }

            this->time_output = 0;
        }
        this->time_output += this->dt;

        // // bi-linear T
        // if (this->flag == 0) {
        //     this->T[0] -= this->v_ch * this->dt;
        // }
        // else {
        //     this->T[0] += this->v_ch * this->dt;
        // }

        // Cosine T
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_ch);

        // check water supply
        x_ub_temp = this->x_ub;
        x_lb_temp = this->x_lb;
        std::sort(x_lb_temp.begin(), x_lb_temp.end());
        std::sort(x_ub_temp.begin(), x_ub_temp.end());

        double d_tot = 0;
        for (int i = 0; i < static_cast<int>(x_lb_temp.size()); i++) {
            d_tot += (x_lb_temp[i] - x_ub_temp[i]);
        }

        for (int j = 0; j < static_cast<int>(x_lb_temp.size()); j++) {
            double d_tot_temp = 0;
            for (int k = j; k < static_cast<int>(x_lb_temp.size()); k++) {
                d_tot_temp += (x_lb_temp[k] - x_ub_temp[k]);
            }
            mw.push_back(d_tot_temp + (this->d - x_ub_temp[j] - d_tot_temp) * this->phi);
            mw.push_back(d_tot_temp + (this->d - x_lb_temp[j] - d_tot_temp) * this->phi - (x_lb_temp[j] - x_ub_temp[j]));
        }

        if (((1 - this->phi) * d_tot < mw.back()) || ((1 - this->phi) * d_tot == mw.back())) {
            this->x_dry = this->d - (1 - this->phi) * d_tot / this->phi;
        }
        else {
            for (int i = 0; i < static_cast<int>(mw.size()) - 1; i++) {
                if ((((1 - this->phi) * d_tot < mw[i]) && ((1 - this->phi) * d_tot > mw[i + 1])) || ((1 - this->phi) * d_tot == mw[i])) {
                    if (i % 2 == 0) {
                        this->x_dry = x_lb_temp[i / 2] - (1 - this->phi) * d_tot + mw[i + 1];
                        break;
                    }
                    else {
                        this->x_dry = x_ub_temp[i / 2 + 1] - ((1 - this->phi) * d_tot - mw[i + 1]) / this->phi;
                        break;
                    }
                }
            }
        }
        mw.clear();

        if ((this->x_dry < this->x0.back()) || (this->x_dry == this->x0.back()) || (this->x_dry < x_lb_temp.back())) {
            this->water_supply = 0;
            this->vl = 0;
            this->solver_ptr->solve(this);
            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }

        //! No ice lens is active
        if (this->ice_active == 1e6) {

            this->vl = 0;
            this->solver_ptr->solve(this);

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }

        // constrain the ice lens boundary
        if (this->Tl[this->ice_active] > this->Tf) {
            if ((this->Tl[this->ice_active] > this->Tm) || (this->Tl[this->ice_active] == this->Tm)) {
                this->vl = 0;
            }
            else {
                this->vl = (this->rho_i * this->L / this->Tm * (this->Tm - this->Tl[this->ice_active])) / (this->mu / this->k0 * (this->x0_active - this->x_lb[this->ice_active]));
            }
            this->solver_ptr->solve(this);

            for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
                if (this->x_lb[i] > this->x_lb[this->ice_active]) {
                    this->x_lb[i] = this->x_lb[i] + this->vl * this->dt;
                    this->x_ub[i] = this->x_ub[i] + this->vl * this->dt;
                    this->idx_lb[i] = this->find_index(this->x_lb[i]);
                    this->idx_ub[i] = this->find_index(this->x_ub[i]);
                }
            }
            this->x_lb[this->ice_active] = this->x_lb[this->ice_active] + this->vl * this->dt;
            this->idx_lb[this->ice_active] = this->find_index(this->x_lb[this->ice_active]);

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }
        
        // calculate heave rate of current ice lens
        this->vl = this->rho_i * this->L / this->Tm * this->FT_int(this->x_lb[this->ice_active]) / this->mu / this->F_mu_int(this->x_lb[this->ice_active]);

        new_ice_lens = false;

        for (int i = 0; i < N_temp; i++) {
            x_temp = this->x_lb[this->ice_active] + (this->x0_active - this->x_lb[this->ice_active]) / (N_temp - 1) * i;
            integral_FT = this->FT_int(x_temp);
            integral_F_mu = this->F_mu_int(x_temp);
            index_x_temp = this->find_index(x_temp);
            if (index_x_temp == this->idx_lb[this->ice_active]) {
                Tx = this->Tl[this->ice_active] + (x_temp - this->x_lb[this->ice_active]) / (this->x[this->idx_lb[this->ice_active] + 1] - this->x_lb[this->ice_active]) * (T_temp[index_x_temp + 1] + 273.15 - this->Tl[this->ice_active]);
            }
            else {
                Tx = this->T[index_x_temp] + (x_temp - this->x[index_x_temp]) / this->dx * (T_temp[index_x_temp + 1] - T_temp[index_x_temp]) + 273.15;
            }

            if (Tx < this->Tf) {
                Fp.push_back(this->mu * this->vl * integral_F_mu + this->c * (1 - this->phi) - this->rho_i * this->L / this->Tm * (integral_FT + (this->Tm - Tx) * (1 - this->phi * this->ice_sat_fun(Tx))));
                x_temp1.push_back(x_temp);
            }
        }

        // get net interparticle force profile above frozen fringe boundary
        if (static_cast<int>(Fp.size()) == 0) {
            this->solver_ptr->solve(this);
            for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
                if (this->x_lb[i] > this->x_lb[this->ice_active]) {
                    this->x_lb[i] = this->x_lb[i] + this->vl * this->dt;
                    this->x_ub[i] = this->x_ub[i] + this->vl * this->dt;
                    this->idx_lb[i] = this->find_index(this->x_lb[i]);
                    this->idx_ub[i] = this->find_index(this->x_ub[i]);
                }
            }
            this->x_lb[this->ice_active] = this->x_lb[this->ice_active] + this->vl * this->dt;
            this->idx_lb[this->ice_active] = this->find_index(this->x_lb[this->ice_active]);

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }

        this->solver_ptr->solve(this);

        // check if new ice lens is initiated
        for (int j = 0; j < static_cast<int>(Fp.size()) - 1; j++) {
            if ((Fp[j] * Fp[j + 1] < 0) || (Fp[j] == 0)) {
                double lens_position = x_temp1[j] + (- Fp[j]) / (Fp[j + 1] - Fp[j]) * (x_temp1[j + 1] - x_temp1[j]);
                new_ice_lens = true; 
                this->x_lb.push_back(lens_position);
                this->x_ub.push_back(lens_position);
                this->idx_lb.push_back(this->find_index(lens_position));
                this->idx_ub.push_back(this->find_index(lens_position));
                this->T_ub.push_back(this->T[this->idx_lb.back()] + (lens_position - this->x[this->idx_lb.back()]) / this->dx * (this->T[this->idx_lb.back() + 1] - this->T[this->idx_lb.back()]) + 273.15);
                this->Tl.push_back(this->T[this->idx_lb.back()] + (lens_position - this->x[this->idx_lb.back()]) / this->dx * (this->T[this->idx_lb.back() + 1] - this->T[this->idx_lb.back()]) + 273.15);
                this->crack_state.push_back(0);
                std::cout << "Ice lens formed at " << lens_position << std::endl;

                break;
            }
        }

        if (new_ice_lens == false) {
            for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
                if (this->x_lb[i] > this->x_lb[this->ice_active]) {
                    this->x_lb[i] = this->x_lb[i] + this->vl * this->dt;
                    this->x_ub[i] = this->x_ub[i] + this->vl * this->dt;
                    this->idx_lb[i] = this->find_index(this->x_lb[i]);
                    this->idx_ub[i] = this->find_index(this->x_ub[i]);
                }
            }
            this->x_lb[this->ice_active] = this->x_lb[this->ice_active] + this->vl * this->dt;
            this->idx_lb[this->ice_active] = this->find_index(this->x_lb[this->ice_active]);

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }
        }
        x_temp1.clear();
        Fp.clear();
    }
}

void Ice_model::kappa_soil() {
    for (int i = 0; i < this->N; i++) {
        this->kappa_s[i] = 5.4e3;
    }
}

double Ice_model::kappa_soil_fun(double x) {
    double k1 = 3.6e3; // thermal conductivity over 0 degree
    double k2 = 5.4e3; // thermal conductivity below -2 degree
    double kappa_sp;
    if (x > 0) {
        kappa_sp = k1;
    }
    else if (x < -2) {
        kappa_sp = k2;
    }
    else {
        kappa_sp = (k2 - k1) * (3 * pow(x, 2) / 4 + pow(x, 3) / 4) + k1;
    }

    return kappa_sp;   
}

double Ice_model::kappa_fun_der(double x) {
    double k1 = 3.6e3;
    double k2 = 5.4e3;
    double kappa_der;
    if (x > 0) {
        kappa_der = 0;
    }
    else if (x < -2) {
        kappa_der = 0;
    }
    else {
        kappa_der = (k2 - k1) * (3 * pow(x, 2) + 6 * x) / 4;
    }

    return kappa_der;
}

double Ice_model::kappa_gradient(double x) {
    double Si = this->ice_sat_fun(x);
    double kappa_g;
    kappa_g = this->kappa_fun_der(x - 273.15) * pow(this->kappa_i, this->phi * Si) * pow(this->kappa_l, this->phi * (1 - Si)) + this->kappa_soil_fun(x - 273.15) * \
     (pow(this->kappa_i, this->phi * Si) * log(this->kappa_i) * this->phi * pow(this->kappa_l, this->phi * (1 - Si)) - pow(this->kappa_i, this->phi * Si) * \
     pow(this->kappa_l, this->phi * (1 - Si)) * log(this->kappa_l) * this->phi) * this->ice_sat_der(x);

     return kappa_g;
}

int Ice_model::find_index(double xx) {
    int index;
    for (int i = 0; i < this->N - 1; i++) {
        if ((xx > this->x[i]) && (xx < this->x[i + 1])) {
            index = i;
            break;
        }
        if (xx == this->x[i]) {
            index = i;
            break;
        }
    }
    return index;
}


void Ice_model::cal_rho_c_eff() {
    std::vector<double> Si = this->ice_sat();
    for (int j = 0; j < this->N; j++) {
        if  (this->label[j] == 2) {
            this->rho_c_eff[j] = this->rho_i * this->cp_i;
        }
        else if (this->label[j] == 3) {
            this->rho_c_eff[j] = this->rho_l * this->cp_l;
        }
        else if (this->label[j] == 1) {
            this->rho_c_eff[j] = Si[j] * this->phi * this->rho_i * this->cp_i + (1 - Si[j]) * this->phi * this->rho_l * this->cp_l + (1 - this->phi) * this->rho_s * this->cp_s;
        }
        else {
            this->rho_c_eff[j] = this->rho_s * this->cp_s * (1 - this->phi);
        }
    }
}

void Ice_model::cal_kappa() {
    std::vector<double> Si = this->ice_sat();
    for (int j = 0; j < this->N; j++) {
        if  (this->label[j] == 2){
            this->kappa[j] = this->kappa_i;
        }
        else if (this->label[j] == 3) {
            this->kappa[j] = this->kappa_l;
        }
        else if (this->label[j] == 1){
            this->kappa[j] = pow(this->kappa_s[j], 1 - this->phi) * pow(this->kappa_i, this->phi * Si[j]) * pow(this->kappa_l, this->phi * (1 - Si[j]));
        }
        else {
            this->kappa[j] = this->kappa_s[j];
        }        
    }
}

void Ice_model::NablaT() {
    for (int i = 0; i < this->N; i++) {
        if (i == 0) {
            this->nablaT[i] = (this->T[1] - this->T[0]) / this->dx;
        }
        else if (i < this->N - 1) {
            this->nablaT[i] = (this->T[i + 1] - this->T[i - 1]) / (2 * this->dx);
        }
        else {
            this->nablaT[i] = (this->T[this->N - 1] - this->T[this->N - 2]) / this->dx;
        }
        for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
            if ((i == this->idx_ub[j]) && (i == 0)) {
                this->nablaT[i] = (this->T_ub[j] - 273.15 - this->T[0]) / this->x_ub[j];
                break;
            }
            else if (i == this->idx_ub[j]) {
                this->nablaT[i] = (this->T_ub[j] - 273.15 - this->T[i - 1]) / (this->dx + this->x_ub[j] - this->x[i]);
                break;
            }
            else if ((i == this->idx_lb[j] + 1) && (i == this->N - 1)) {
                this->nablaT[i] = (this->T[this->N - 1] + 273.15 - this->Tl[j]) / (this->x[i] - this->x_lb[j]);
                break;
            }
            else if (i == this->idx_lb[j] + 1) {
                this->nablaT[i] = (this->T[i + 1] + 273.15 - this->Tl[j]) / (this->x[i] - this->x_lb[j] + this->dx);
                break;
            }
        }
    }
}

std::vector<double> Ice_model::ice_sat() {
    std::vector<double> si_temp(this->N, 0);
    for (int i = 0; i < this->N; i++) {
        si_temp[i] = this->ice_sat_fun(this->T[i] + 273.15);
    }

    return si_temp;
}

double Ice_model::ice_sat_fun(double x) {
    if (x < this->Tf) {
        return 1 - pow((this->Tm - this->Tf) / (this->Tm - x), 2);
    }
    // else if (x < this->Tf + 1e-3) {
    //     double si = 1 - pow((this->Tm - this->Tf) / (this->Tm - (this->Tf - 1e-3)), 2);
    //     double R = (pow(si, 2) + pow(2e-3, 2)) / (2 * si);
    //     return R - sqrt(pow(R, 2) - pow(x-(this->Tf + 1e-3), 2));
    // }
    else {
        return 0;
    }
}

double Ice_model::ice_sat_der(double x) {
    if (x < this->Tf) {
        return 2 * pow(this->Tm - this->Tf, 2) * pow(x - this->Tm, -3);
    }
    // else if (x < this->Tf + 1e-3) {
    //     double si = 1 - pow((this->Tm - this->Tf) / (this->Tm - (this->Tf - 1e-3)), 2);
    //     double R = (pow(si, 2) + pow(2e-3, 2)) / (2 * si);
    //     return (x - (this->Tf + 1e-3)) / sqrt(pow(R, 2) - pow(x - (this->Tf + 1e-3), 2));
    // }
    else {
        return 0;
    }
}

double Ice_model::FT_int(double x) {
    // Calculate integral occurring in expression of FT
    double integral = 0;
    int N_temp = 20;
    std::vector<double> x_temp(N_temp, 0);
    for (int i = 0; i < N_temp; i++) {
        x_temp[i] = x + i * (this->x0_active - x) / (N_temp - 1);
    }
    // int index = this->find_index(x);
    for (int i = 0; i < N_temp - 1; i++) {
        integral += (this->FT_int_fun(x_temp[i]) + this->FT_int_fun(x_temp[i + 1])) / 2 * (this->x0_active - x) / (N_temp - 1);
    }

    return integral;
}

double Ice_model::FT_int_fun(double x) {
    double Tx; // temperature at x
    double nablaTx;
    int index = this->find_index(x);
    if (this->idx_lb.empty()) {
        Tx = this->T[index] + (x - index * this->dx) / this->dx * (this->T[index + 1] - this->T[index]) + 273.15;
        nablaTx = this->nablaT[index] + (this->nablaT[index + 1] - this->nablaT[index]) * (x - index * this->dx) / this->dx;
    }
    else {
        double nablaTl = (this->T[this->idx_lb[this->ice_active] + 1] + 273.15 - this->Tl[this->ice_active]) / (this->x[index + 1] - this->x_lb[this->ice_active]);
        if (index == this->idx_lb[this->ice_active]) {
            Tx = this->Tl[this->ice_active] + (x - this->x_lb[this->ice_active]) / (this->x[index + 1] - this->x_lb[this->ice_active]) * (this->T[index + 1] + 273.15 - this->Tl[this->ice_active]);
            nablaTx = nablaTl + (this->nablaT[index + 1] - nablaTl) * (x - this->x_lb[this->ice_active]) / (this->x[index + 1] - this->x_lb[this->ice_active]);
        }
        else {
            Tx = this->T[index] + (x - index * this->dx) / this->dx * (this->T[index + 1] - this->T[index]) + 273.15;
            nablaTx = this->nablaT[index] + (this->nablaT[index + 1] - this->nablaT[index]) * (x - index * this->dx) / this->dx;
        }
    }

    return nablaTx * (1 - this->phi * this->ice_sat_fun(Tx));
}

double Ice_model::F_mu_int(double x) {
    //Calculate the integral occurring in expression of F_mu
    // int index = this->find_index(x);
    double integral = 0;
    int N_temp = 20;
    std::vector<double> x_temp(N_temp, 0);
    for (int i = 0; i < N_temp; i++) {
        x_temp[i] = x + i * (this->x0_active - x) / (N_temp - 1);
    }
    
    for (int i = 0; i < N_temp - 1; i++) {
        integral += (this->F_mu_int_fun(x_temp[i]) + this->F_mu_int_fun(x_temp[i + 1])) / 2 * (this->x0_active - x) / (N_temp - 1);
    }

    return integral;

}

double Ice_model::F_mu_int_fun(double x) {
    int index = this->find_index(x);
    double Tx;
    Tx = this->T[index] + (x - index * this->dx) / this->dx * (this->T[index + 1] - this->T[index]) + 273.15;   
    double Si_x = this->ice_sat_fun(Tx);
    double kx = this->k0 * pow(1 - Si_x, 2);
    return pow(1 - this->phi * Si_x, 2) / kx;
}

int Ice_model::ice_lens_active() {
    int idx = 0;
    bool ice_lens_active = true;
    bool full_ice = false;

    for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
        if ((this->crack_state[i] == 0) || ((this->crack_state[i] == 1) && (this->Tl[i] < this->Tm))) {
            idx = i;
            full_ice = true;
        }
    }

    for (int j = idx + 1; j < static_cast<int>(this->x_lb.size()); j++) {
        if (this->crack_state[j] == 1) {
            ice_lens_active = false;
            break;
        }
    }

    if (((idx != static_cast<int>(this->x_lb.size()) - 1) && (ice_lens_active == false)) || (full_ice == false)) {
        idx = 1e6;
    }

    return idx;
}
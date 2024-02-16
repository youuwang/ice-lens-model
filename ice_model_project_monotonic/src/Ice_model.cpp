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
        sstr >> this->v_cd;
    }
  }
}


void Ice_model::solve(ofstream& output_T, ofstream& output_ice_lens) {

    ofstream time_N_lens(this->out_dir  + "/t_N.dat");
    ofstream ub_lb(this->out_dir  + "/ub_lb_t.dat");
    ofstream time_xf2(this->out_dir  + "/t_xf2.dat");
    ofstream time_xf3(this->out_dir  + "/t_xf3.dat");
    ofstream Tu_Tl(this->out_dir  + "/Tu_Tl.dat");
    ofstream time1(this->out_dir  + "/t1.dat");
    ofstream time_dry(this->out_dir  + "/t_dry.dat");
    this->solve1(output_T);
    std::cout << "Phase 1 finished ..." << std::endl;
    if (this->T[0] < this->T1) {
        std::cout << "No frozen fringe formed " << std::endl;
        return;
    }

    this->solve2(output_T);
    std::cout << "Phase 2 finished ..." << std::endl;
    if (this->T[0] < this->T1) {
        std::cout << "Frozen fringe involved but no ice lens formed " << this->v << std::endl;
        return;
    }

    this->solve3(output_T);
    std::cout << "Phase 3 finished and ice lenses are:  " << std::endl;

    for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
        std::cout << this->x_ub[i] << "\t" << this->x_lb[i] << std::endl;
        output_ice_lens << this->x_ub[i] << "\t" << this->x_lb[i] << endl;
    }

    for (int m = 0; m < static_cast<int>(this->xf_t2.size()); m++) {
        time_xf2 << this->time_t2[m] << "\t" << this->xf_t2[m] << endl;
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

    for (int r = 0; r < static_cast<int>(this->x_dry_t.size()); r++) {
        time_dry << this->time_t3[r] << "\t" << this->x_dry_t[r] << endl;
    }

    time_N_lens.close();
    ub_lb.close();
    time_xf2.close();
    time_xf3.close();
    Tu_Tl.close();
    time1.close();
    time_dry.close();

}


void Ice_model::solve1(ofstream& output_T) {

    this->solver_ptr = std::make_shared<FDM_solver>(this->d, this->N);

    // Solve the heat equation using the finite difference method
    while (true) {
        this->time += this->dt;

        if (this->T[0] < (this->Tf - 273.15)) {
            return;
        }

        if (this->time > this->v_cd / 2) {
            return;
        }

        //output results
        if (fabs(this->time_output - 0.01) < 1e-8) {
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
        
        // // linear T
        // this->T[0] -= this->v_cd * this->dt; 

        // cosine T
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_cd);

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

    this->solver_ptr = std::make_shared<FDM_ff>(this->d, this->N);

    ofstream F_p(this->out_dir  + "/F_p.dat");

    // Solve the heat equation using the finite difference method
    while (true) {
        this->time += this->dt;

        // Criterion for ending the second phase (ice lens model comes into play)
        if (this->time > this->v_cd / 2) {
            return;
        }

        // locate Tm
        for (int i = 0; i < this->N - 1; i++) {
            if ((this->T[i] * this->T[i + 1] < 0) || (this->T[i] == 0)) {
                this->x0 = this->x[i] + this->dx * (-this->T[i]) / (this->T[i + 1] - this->T[i]);
                break;
            }
        }

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
        if (fabs(this->time_output - 0.01) < 1e-8) {
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

        // // linear T
        // this->T[0] -= this->v_cd * this->dt;
        // // Criterion for ending the second phase (ice lens model comes into play)
        // if (this->T[0] < this->T1) {
        //     return;
        // }        

        // cosine T
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_cd);

        // call the solver
        this->solver_ptr->solve(this);

        for (int j = 0; j < N_temp - 1; j++) {
            if ((Fp[j] * Fp[j + 1] < 0) || (Fp[j] == 0)) {
                this->lens_ini = x_temp[j] + (-Fp[j]) / (Fp[j + 1] - Fp[j]) * this->xf / (N_temp - 1);
                std::cout << "Ice lens formed at " << this->lens_ini << " at time " << this->time << std::endl;
                this->x_lb.push_back(this->lens_ini);
                this->x_ub.push_back(this->lens_ini);
                this->idx_lb.push_back(this->find_index(this->lens_ini));
                this->idx_ub.push_back(this->find_index(this->lens_ini));

                int ice_active_ini = this->find_index(this->lens_ini);
                this->Tl.push_back(this->T[ice_active_ini] + (this->lens_ini - this->x[ice_active_ini]) / this->dx * (this->T[ice_active_ini + 1] - this->T[ice_active_ini]) + 273.15);
                this->T_ub.push_back(this->T[ice_active_ini] + (this->lens_ini - this->x[ice_active_ini]) / this->dx * (this->T[ice_active_ini + 1] - this->T[ice_active_ini]) + 273.15);

                this->vl = this->rho_i * this->L / this->Tm * this->FT_int(this->x_lb.back()) / this->mu / this->F_mu_int(this->x_lb.back());
                std::cout << "Heave rate is " << this->vl << std::endl;
                return;
            }
        }
    }
    F_p.close();
}

void Ice_model::solve3(ofstream& output_T) {
    double integral_FT; // integral inside the expression of FT
    double integral_F_mu; // integral inside the expression of F_mu
    int N_temp = 20; // size of x_temp
    std::vector<double> Fp(N_temp, 0); // net interparticle force at each node of interest
    std::vector<double> x_temp(N_temp, 0); // temporary x for calculation of Fp
    int index_x_temp; // index of element where x_temp is located
    double Tx; // temperature at x_temp
    this->solver_ptr = std::make_shared<FDM_icelens>(this->d, this->N);
    std::vector<double> T_temp;
    bool new_ice_lens;

    while (true) {
        this->time += this->dt;

        // stopping criterion
        if (this->time > this->v_cd / 2) {
            return;
        }

        // check water supply
        double d_tot = 0;
        for (int i = 0; i < static_cast<int>(this->x_lb.size()); i++) {
            d_tot += (this->x_lb[i] - this->x_ub[i]);
        }
        this->x_dry = this->d - d_tot * (1 - this->phi) / this->phi;

        // label nodes
        for (int i = 0; i < this->N; i++) {
            if ((this->x[i] > this->x_dry) || (this->x[i] == this->x_dry)) {
                this->label[i] = 3;
            }
            else {
                for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                    if ((this->x[i] > this->x_ub[j]) && (this->x[i] < this->x_lb[j])) {
                        this->label[i] = 2;
                        break;
                    }
                }
            }
        }

        // locate Tf
        for (int i = 0; i < this->N - 1; i++) {
            if ((((this->Tf - 273.15) - this->T[i]) * ((this->Tf - 273.15) - this->T[i + 1]) < 0) || ((this->Tf - 273.15) - this->T[i] == 0)) {
                this->xf = this->x[i] + this->dx * (this->Tf - 273.15 - this->T[i]) / (this->T[i + 1] - this->T[i]);
                break;
            }    
        }

        // locate Tm
        for (int i = 0; i < this->N - 1; i++) {
            if ((this->T[i] * this->T[i + 1] < 0) || (this->T[i] == 0)) {
                this->x0 = this->x[i] + this->dx * (- this->T[i]) / (this->T[i + 1] - this->T[i]);
                break;
            }
        }

        if (this->x_lb.back() > this->d) {
            std::cout << "impossible over d" << std::endl;
            return;
        }

        T_temp = this->T;

        //! output results
        if (fabs(this->time_output - 0.01) < 1e-8) {
            for (int i = 0; i < this->N; i++) {
                if (i < this->N - 1) {
                    output_T << T_temp[i] << "\t";
                }
                else {
                    output_T << T_temp[i] << endl;
                }
            }

            for (int i = 0; i < static_cast<int>(this->Tl.size()); i++) {
                std::cout << this->T_ub[i] << "\t" << this->Tl[i] <<  "\t" << this->vl << std::endl;
            }

            this->time_t3.push_back(this->time);
            this->xf_t3.push_back(this->xf);
            this->x_dry_t.push_back(this->x_dry);
            for (int j = 0; j < static_cast<int>(this->x_ub.size()); j++) {
                this->x_ub_t.push_back(this->x_ub[j]);
                this->x_lb_t.push_back(this->x_lb[j]);
                this->T_ub_t.push_back(this->T_ub[j]);
                this->Tl_t.push_back(this->Tl[j]);
            }
            this->N_lens.push_back(static_cast<int>(this->x_ub.size()));
            this->time_output = 0;
        }
        this->time_output += this->dt;

        // // linear T
        // this->T[0] -= this->v_cd * this->dt;

        // // stopping criterion
        // if (this->T[0] < this->T1) {
        //     return;
        // }

        // cosine T
        this->T[0] = (this->T1 + this->T2) / 2 + (this->T2 - this->T1) / 2 * cos(2 * M_PI * this->time / this->v_cd);

        if ((this->x_dry < this->x0) || (this->x_dry == this->x0)) {
            std::cout << "No water supply" << std::endl;
            this->vl = 0;
            this->solver_ptr->solve(this);
            
            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }
        

        if (this->x_lb.back() > this->xf) {

            if ((this->x_lb.back() > this->x0) || (this->x_lb.back() == this->x0)) {
                this->vl = 0;
            }
            else {
                this->vl = (this->rho_i * this->L / this->Tm * (this->Tm - this->Tl.back())) / (this->mu / this->k0 * (this->x0 - this->x_lb.back()));
            }

            this->solver_ptr->solve(this);
            this->x_lb.back() = this->x_lb.back() + this->vl * this->dt;
            this->idx_lb.back() = this->find_index(this->x_lb.back());

            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }

            continue;
        }

        // calculate growth rate of current ice lens
        this->vl = this->rho_i * this->L / this->Tm * this->FT_int(this->x_lb.back()) / this->mu / this->F_mu_int(this->x_lb.back());

        if (this->vl < 0) {
            std::cout << "impossible" << std::endl;
            this->vl = 0;
        }

        // get net interparticle force profile above frozen fringe boundary
        for (int i = 0; i < N_temp; i++) {
            x_temp[i] = this->x_lb.back() + (this->xf - this->x_lb.back()) / (N_temp - 1) * i;
        }
        for (int j = 0; j < N_temp; j++) {
            integral_FT = this->FT_int(x_temp[j]);
            integral_F_mu = this->F_mu_int(x_temp[j]);
            index_x_temp = this->find_index(x_temp[j]);
            if (index_x_temp == this->idx_lb.back()) {
                Tx = this->Tl.back() + (x_temp[j] - this->x_lb.back()) / (this->x[this->idx_lb.back() + 1] - this->x_lb.back()) * (this->T[index_x_temp + 1] + 273.15 - this->Tl.back());
            }
            else {
                Tx = this->T[index_x_temp] + (x_temp[j] - this->x[index_x_temp]) / this->dx * (this->T[index_x_temp + 1] - this->T[index_x_temp]) + 273.15;
            }
            Fp[j] = this->mu * this->vl * integral_F_mu + this->c * (1 - this->phi) - this->rho_i * this->L / this->Tm * (integral_FT + (this->Tm - Tx) * (1 - this->phi * this->ice_sat_fun(Tx)));
        }

        // get idx during iterations
        this->idx = this->find_index(this->xf);

        // update temperature profile
        this->solver_ptr->solve(this);

        new_ice_lens = false;

        // check if new ice lens is initiated
        for (int j = 0; j < N_temp - 1; j++) {
            if ((Fp[j] * Fp[j + 1] < 0) || (Fp[j] == 0)) {
                double lens_position = x_temp[j] + (- Fp[j]) / (Fp[j + 1] - Fp[j]) * (this->xf - this->x_lb.back()) / (N_temp - 1);
                new_ice_lens = true;
                this->x_lb.push_back(lens_position);
                this->x_ub.push_back(lens_position);
                this->idx_lb.push_back(this->find_index(lens_position));
                this->idx_ub.push_back(this->find_index(lens_position));
                this->Tl.push_back(this->T[this->idx_lb.back()] + (this->T[this->idx_lb.back() + 1] - this->T[this->idx_lb.back()]) / this->dx * (lens_position - this->x[this->idx_lb.back()]) + 273.15);
                this->T_ub.push_back(this->T[this->idx_lb.back()] + (this->T[this->idx_lb.back() + 1] - this->T[this->idx_lb.back()]) / this->dx * (lens_position - this->x[this->idx_lb.back()]) + 273.15);
                std::cout << "Ice lens formed at " << lens_position << "\t" << this->vl << "\t" << this->Tl.back() << "\t" << this->T_ub.back() << std::endl; 

                break;
            }
        }

        //! update ice lens boundary
        if (new_ice_lens == false) {
            this->x_lb.back() = this->x_lb.back() + this->vl * this->dt;
            this->idx_lb.back() = this->find_index(this->x_lb.back());
            for (int j = 0; j < static_cast<int>(this->x_lb.size()); j++) {
                this->T_ub[j] = 273.15 + this->T[this->idx_ub[j]] + (this->x_ub[j] - this->x[this->idx_ub[j]]) / this->dx * (this->T[this->idx_ub[j] + 1] - this->T[this->idx_ub[j]]);
                this->Tl[j] = 273.15 + this->T[this->idx_lb[j]] + (this->x_lb[j] - this->x[this->idx_lb[j]]) / this->dx * (this->T[this->idx_lb[j] + 1] - this->T[this->idx_lb[j]]);
            }
        }
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
            this->rho_c_eff[j] = this->rho_s * this->cp_s;
        }
        else {
            this->rho_c_eff[j] =  Si[j] * this->phi * this->rho_i * this->cp_i + (1 - Si[j]) * this->phi * this->rho_l * this->cp_l + (1 - this->phi) * this->rho_s * this->cp_s;
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
            this->kappa[j] = this->kappa_s[j];
        }
        else {
            // this->kappa[j] = pow(this->kappa_s[j], 1 - this->phi) * pow(this->kappa_l, this->phi);
            this->kappa[j] = pow(this->kappa_s[j], 1 - this->phi) * pow(this->kappa_i, this->phi * Si[j]) * pow(this->kappa_l, this->phi * (1 - Si[j]));
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

    else {
        return 0;
    }
}

double Ice_model::ice_sat_der(double x) {
    if (x < this->Tf) {
        return 2 * pow(this->Tm - this->Tf, 2) * pow(x - this->Tm, -3);
    }

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
        x_temp[i] = x + i * (this->xf - x) / (N_temp - 1);
    }
    // int index = this->find_index(x);
    integral += (this->Tm - this->Tf);
    for (int i = 0; i < N_temp - 1; i++) {
        integral += (this->FT_int_fun(x_temp[i]) + this->FT_int_fun(x_temp[i + 1])) / 2 * (this->xf - x) / (N_temp - 1);
    }

    return integral;
}

double Ice_model::FT_int_fun(double xx) {
    double Tx; // temperature at x
    double nablaTx;
    int index = this->find_index(xx);
    if (this->idx_lb.empty()) {
        Tx = this->T[index] + (xx - index * this->dx) / this->dx * (this->T[index + 1] - this->T[index]) + 273.15;
        nablaTx = this->nablaT[index] + (this->nablaT[index + 1] - this->nablaT[index]) * (xx - index * this->dx) / this->dx;
    }
    else {
        double nablaTl = (this->T[this->idx_lb.back() + 1] + 273.15 - this->Tl.back()) / (this->x[index + 1] - this->x_lb.back());
        if (index == this->idx_lb.back()) {
            Tx = this->Tl.back() + (xx - this->x_lb.back()) / (this->x[index + 1] - this->x_lb.back()) * (this->T[index + 1] + 273.15 - this->Tl.back());
            nablaTx = nablaTl + (this->nablaT[index + 1] - nablaTl) * (xx - this->x_lb.back()) / (this->x[index + 1] - this->x_lb.back());
        }
        else {
            Tx = this->T[index] + (xx - index * this->dx) / this->dx * (this->T[index + 1] - this->T[index]) + 273.15;
            nablaTx = this->nablaT[index] + (this->nablaT[index + 1] - this->nablaT[index]) * (xx - index * this->dx) / this->dx;
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
        x_temp[i] = x + i * (this->xf - x) / (N_temp - 1);
    }

    integral += (this->x0 - this->xf) / this->k0;
    
    for (int i = 0; i < N_temp - 1; i++) {
        integral += (this->F_mu_int_fun(x_temp[i]) + this->F_mu_int_fun(x_temp[i + 1])) / 2 * (this->xf - x) / (N_temp - 1);
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
#include "FDM_icelens.hpp"

void FDM_icelens::solve(Ice_model* ptr) {

    // calculate rho_c_eff for each node
    ptr->cal_rho_c_eff();

    // calculate kappa for each node
    ptr->kappa_soil();
    ptr->cal_kappa();

    std::vector<double> T_temp = ptr->T;

    if (ptr->ice_active != 1e6) {
        double xl = ptr->x_lb[ptr->ice_active];
        int idx_xl = ptr->idx_lb[ptr->ice_active];
        double si_active = ptr->ice_sat_fun(ptr->Tl[ptr->ice_active]);

        for (int i = 1; i < this->N - 1; i++) {
            if (i == idx_xl) {
                ptr->T[i] = T_temp[i] + this->dt * (((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i - 1] + ptr->kappa[i]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) + \
                ptr->rho_i * ptr->L * ptr->vl * (1 - ptr->phi * si_active) * (ptr->x[idx_xl + 1] - xl) / ptr->dx) / ptr->rho_c_eff[i];
            }
            else if (i == idx_xl + 1) {
                ptr->T[i] = T_temp[i] + this->dt * (((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i - 1] + ptr->kappa[i]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) + \
                ptr->rho_i * ptr->L * ptr->vl * (1 - ptr->phi * si_active) * (xl - ptr->x[idx_xl]) / ptr->dx) / ptr->rho_c_eff[i];        
            }
            else {
                ptr->T[i] = T_temp[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i - 1] + ptr->kappa[i]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) / ptr->rho_c_eff[i];
            }
        }      
    }
    else {
        for (int i = 1; i < this->N - 1; i++) {
            ptr->T[i] = T_temp[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i - 1] + ptr->kappa[i]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) / ptr->rho_c_eff[i];
        }
    }

/*
    double Tl_active;

    if (ptr->ice_active != 1e6) {
        double xl = ptr->x_lb[ptr->ice_active];
        double Tl_active_next = 0;
        int idx_xl = ptr->find_index(xl);
        Tl_active = 273.15 + ptr->T[idx_xl] + (ptr->T[idx_xl + 1] - ptr->T[idx_xl]) * (xl - ptr->x[idx_xl]) / this->dx;
        int idx_loc;
        int n_loc = 1000;
        double dx_loc = this->dx / (n_loc - 1);
        std::vector<double> x_loc(n_loc, 0);
        std::vector<double> T_loc(n_loc, 0);
        double kappa_loc;
        double kappa_s_loc;
        double si_loc;
        double f_value_loc;
        double f_der_loc;
        for (int i = 0; i < n_loc; i++) {
            x_loc[i] = ptr->x[idx_xl] + dx_loc * i;
            T_loc[i] = ptr->T[idx_xl] + (ptr->T[idx_xl + 1] - ptr->T[idx_xl]) * dx_loc * i / this->dx;
        }

        for (int j = 0; j < n_loc - 1; j++) {
            if (((xl > x_loc[j]) && (xl < x_loc[j + 1])) || (xl == x_loc[j])) {
                idx_loc = j;
                break;
            }
        }

        double dx_l_loc = xl - x_loc[idx_loc];
        double dx_r_loc = x_loc[idx_loc + 1] - xl;
        while (fabs(Tl_active - Tl_active_next) > 1e-3) {
            Tl_active_next = Tl_active;
            si_loc = ptr->ice_sat_fun(Tl_active);
            kappa_s_loc = ptr->kappa_soil_fun(Tl_active - 273.15);
            kappa_loc = pow(kappa_s_loc, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_loc) * pow(ptr->kappa_l, ptr->phi * (1 - si_loc));

            f_value_loc = ptr->kappa_i * (Tl_active - (T_loc[idx_loc] + 273.15)) / dx_l_loc - kappa_loc * (T_loc[idx_loc + 1] + 273.15 - Tl_active) / dx_r_loc - ptr->rho_i * ptr->L * ptr->vl * (1 - ptr->phi * si_loc);

            f_der_loc = ptr->kappa_i / dx_l_loc + kappa_loc / dx_r_loc - ptr->kappa_gradient(Tl_active) * (T_loc[idx_loc + 1] + 273.15 - Tl_active) / dx_r_loc + ptr->rho_i * ptr->L * ptr->vl * ptr->phi * ptr->ice_sat_der(Tl_active);

            Tl_active = Tl_active - f_value_loc / f_der_loc;
        }
    }

    for (int i = 0; i < static_cast<int>(ptr->x_lb.size()); i++) {
        double dx_l_lb = ptr->x_lb[i] - ptr->x[ptr->idx_lb[i]];
        double dx_r_lb = this->dx - dx_l_lb;
        double dx_l_ub = ptr->x_ub[i] - ptr->x[ptr->idx_ub[i]];
        double dx_r_ub = this->dx - dx_l_ub;
        double delta = ptr->x_lb[i] - ptr->x_ub[i];

        double si_r_lb;
        double Tl_temp = (ptr->T[ptr->idx_lb[i]] + ptr->T[ptr->idx_lb[i] + 1]) / 2 + 273.15;
        double Tl_temp_next = 0;
        double T_ub_temp = (ptr->T[ptr->idx_ub[i]] + ptr->T[ptr->idx_ub[i] + 1]) / 2 + 273.15;
        double T_ub_temp_next = 0;
        double kappa_s_r_lb;
        double kappa_r_lb;
        double kappa_s_l_ub;
        double kappa_l_ub;
        double si_l_ub;
        double g_value;
        double g_der;
        double f_value;
        double f_der;
        int inice_0 = 0;

        for (int j = 0; j < static_cast<int>(ptr->x0.size()); j++) {
            if ((ptr->x_lb[i] > ptr->x0[j]) && (ptr->x_ub[i] < ptr->x0[j])) {
                inice_0 = 1;
                break;
            }
        }

        if (inice_0 == 1) {
            if (ptr->T_ub[i] - 273.15 < 0) {
                while (fabs(T_ub_temp - T_ub_temp_next) > 1e-3) {
                    T_ub_temp_next = T_ub_temp;
                    kappa_s_l_ub = ptr->kappa_soil_fun(T_ub_temp - 273.15);
                    si_l_ub = ptr->ice_sat_fun(T_ub_temp);
                    kappa_l_ub = pow(kappa_s_l_ub, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_l_ub) * pow(ptr->kappa_l, ptr->phi * (1 - si_l_ub));

                    g_value = kappa_l_ub * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub - ptr->kappa_i * (ptr->T[ptr->idx_ub[i] + 1] + 273.15 - T_ub_temp) / dx_r_ub;

                    g_der = ptr->kappa_gradient(T_ub_temp) * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub + kappa_l_ub / dx_l_ub + ptr->kappa_i / dx_r_ub;

                    T_ub_temp = T_ub_temp - g_value / g_der;
                }
                if (i != ptr->ice_active) {
                    while (fabs(Tl_temp - Tl_temp_next) > 1e-3) {
                        Tl_temp_next = Tl_temp;
                        si_r_lb = ptr->ice_sat_fun(Tl_temp);
                        kappa_s_r_lb = ptr->kappa_soil_fun(Tl_temp - 273.15);
                        kappa_r_lb = pow(kappa_s_r_lb, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_r_lb) * pow(ptr->kappa_l, ptr->phi * (1 - si_r_lb));

                        f_value = kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - ptr->kappa_l * (Tl_temp - (ptr->T[ptr->idx_lb[i]] + 273.15)) / dx_l_lb;

                        f_der = ptr->kappa_gradient(Tl_temp) * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_r_lb / dx_r_lb - ptr->kappa_l / dx_l_lb;

                        Tl_temp = Tl_temp - f_value / f_der;
                    }
                }                
            }
            
            else {
                while (fabs(T_ub_temp - T_ub_temp_next) > 1e-3) {
                    T_ub_temp_next = T_ub_temp;
                    kappa_s_l_ub = ptr->kappa_soil_fun(T_ub_temp - 273.15);
                    si_l_ub = ptr->ice_sat_fun(T_ub_temp);
                    kappa_l_ub = pow(kappa_s_l_ub, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_l_ub) * pow(ptr->kappa_l, ptr->phi * (1 - si_l_ub));

                    g_value = kappa_l_ub * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub - ptr->kappa_l * (ptr->T[ptr->idx_ub[i] + 1] + 273.15 - T_ub_temp) / dx_r_ub;

                    g_der = ptr->kappa_gradient(T_ub_temp) * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub + kappa_l_ub / dx_l_ub + ptr->kappa_l / dx_r_ub;

                    T_ub_temp = T_ub_temp - g_value / g_der;
                }
                if (i != ptr->ice_active) {
                    while (fabs(Tl_temp - Tl_temp_next) > 1e-3) {
                        Tl_temp_next = Tl_temp;
                        si_r_lb = ptr->ice_sat_fun(Tl_temp);
                        kappa_s_r_lb = ptr->kappa_soil_fun(Tl_temp - 273.15);
                        kappa_r_lb = pow(kappa_s_r_lb, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_r_lb) * pow(ptr->kappa_l, ptr->phi * (1 - si_r_lb));

                        f_value = kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - ptr->kappa_i * (Tl_temp - (ptr->T[ptr->idx_lb[i]] + 273.15)) / dx_l_lb;

                        f_der = ptr->kappa_gradient(Tl_temp) * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_r_lb / dx_r_lb - ptr->kappa_i / dx_l_lb;

                        Tl_temp = Tl_temp - f_value / f_der;
                    }
                }
            }
        }

        else {
            double kappa_crack;
            if (ptr->crack_state[i] == 0) {
                kappa_crack = ptr->kappa_i;
            }
            else {
                kappa_crack = ptr->kappa_l;
            }
            if (ptr->idx_lb[i] == ptr->idx_ub[i]) {
                while (fabs(Tl_temp - Tl_temp_next) > 1e-3) {
                    Tl_temp_next = Tl_temp;
                    T_ub_temp = Tl_temp - kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb * delta / kappa_crack;
                    kappa_s_l_ub = ptr->kappa_soil_fun(T_ub_temp - 273.15);
                    si_l_ub = ptr->ice_sat_fun(T_ub_temp);
                    kappa_l_ub = pow(kappa_s_l_ub, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_l_ub) * pow(ptr->kappa_l, ptr->phi * (1 - si_l_ub));

                    g_value = kappa_l_ub * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub;
                    
                    g_der = (ptr->kappa_gradient(T_ub_temp) * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub + kappa_l_ub / dx_l_ub) * (1 - delta / kappa_crack * (ptr->kappa_gradient(Tl_temp) * \
                    (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_r_lb / dx_r_lb));
                    
                    f_der = g_der - (ptr->kappa_gradient(Tl_temp) * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_r_lb / dx_r_lb);
                    
                    f_value = g_value - kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb;

                    Tl_temp = Tl_temp - f_value / f_der;
                }
                T_ub_temp = Tl_temp - kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb * delta / kappa_crack;             
            }
            else {
                while (fabs(T_ub_temp - T_ub_temp_next) > 1e-3) {
                    T_ub_temp_next = T_ub_temp;
                    kappa_s_l_ub = ptr->kappa_soil_fun(T_ub_temp - 273.15);
                    si_l_ub = ptr->ice_sat_fun(T_ub_temp);
                    kappa_l_ub = pow(kappa_s_l_ub, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_l_ub) * pow(ptr->kappa_l, ptr->phi * (1 - si_l_ub));

                    g_value = kappa_l_ub * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub - kappa_crack * (ptr->T[ptr->idx_ub[i] + 1] + 273.15 - T_ub_temp) / dx_r_ub;

                    g_der = ptr->kappa_gradient(T_ub_temp) * (T_ub_temp - ptr->T[ptr->idx_ub[i]] - 273.15) / dx_l_ub + kappa_l_ub / dx_l_ub + kappa_crack / dx_r_ub;

                    T_ub_temp = T_ub_temp - g_value / g_der;
                }
                
                if (i != ptr->ice_active) {
                    while (fabs(Tl_temp - Tl_temp_next) > 1e-3) {
                        Tl_temp_next = Tl_temp;
                        si_r_lb = ptr->ice_sat_fun(Tl_temp);
                        kappa_s_r_lb = ptr->kappa_soil_fun(Tl_temp - 273.15);
                        kappa_r_lb = pow(kappa_s_r_lb, 1 - ptr->phi) * pow(ptr->kappa_i, ptr->phi * si_r_lb) * pow(ptr->kappa_l, ptr->phi * (1 - si_r_lb));

                        f_value = kappa_r_lb * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_crack * (Tl_temp - (ptr->T[ptr->idx_lb[i]] + 273.15)) / dx_l_lb;

                        f_der = ptr->kappa_gradient(Tl_temp) * (ptr->T[ptr->idx_lb[i] + 1] + 273.15 - Tl_temp) / dx_r_lb - kappa_r_lb / dx_r_lb - kappa_crack / dx_l_lb;

                        Tl_temp = Tl_temp - f_value / f_der;
                    }
                }
            }
        }

        ptr->T_ub[i] = T_ub_temp; 

        ptr->Tl[i] = Tl_temp;

    }

    if (ptr->ice_active != 1e6) {
        ptr->Tl[ptr->ice_active] = Tl_active;
    }
*/

    // calculate temperature gradient at each node
    ptr->NablaT();
}
#include "FDM_ff.hpp"

void FDM_ff::solve(Ice_model* ptr) {
    
    // Prepare rho_c_eff for each element
    ptr->cal_rho_c_eff();

    // Prepare thermal conductivity for each node
    ptr->cal_kappa();

    std::vector<double> T_temp = ptr->T;

    for (int i = 1; i < this->N - 1; i++) {
        ptr->T[i] = T_temp[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) / ptr->dx[i + 1] - (ptr->kappa[i] + ptr->kappa[i - 1]) / 2 * (T_temp[i] - T_temp[i - 1]) / ptr->dx[i]) / ((ptr->dx[i] + ptr->dx[i + 1]) / 2) / ptr->rho_c_eff[i];
    }

    // calculate temperature gradients at each node
    ptr->NablaT();
}
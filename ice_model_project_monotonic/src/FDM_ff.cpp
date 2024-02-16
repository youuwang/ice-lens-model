#include "FDM_ff.hpp"

void FDM_ff::solve(Ice_model* ptr) {
    
    // Prepare rho_c_eff for each element
    ptr->cal_rho_c_eff();

    // Prepare thermal conductivity for each node
    ptr->kappa_soil();
    ptr->cal_kappa();

    std::vector<double> T_temp = ptr->T;

    for (int i = 1; i < this->N - 1; i++) {
        ptr->T[i] = T_temp[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i] + ptr->kappa[i - 1]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) / ptr->rho_c_eff[i];
    }

    // locate Tf
    for (int i = 0; i < this->N - 1; i++) {
        if (((ptr->Tf - 273.15) > ptr->T[i]) && ((ptr->Tf - 273.15) < ptr->T[i + 1])) {
            ptr->xf = ptr->x[i] + (ptr->x[i + 1] - ptr->x[i]) * ((ptr->Tf - 273.15) - ptr->T[i]) / (ptr->T[i + 1] - ptr->T[i]);
            break;
        }
    }

    // calculate temperature gradients at each node
    ptr->NablaT();

    // get idx during iterations
    ptr->idx = ptr->find_index(ptr->xf);
}
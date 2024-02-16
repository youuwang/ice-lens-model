#include "FDM_solver.hpp"

void FDM_solver::solve(Ice_model* ptr) {

    // calculate temperature gradient at each node
    ptr->NablaT();

    // Prepare rho_c_eff for each element
    ptr->cal_rho_c_eff();

    // Prepare thermal conductivity for each node
    ptr->kappa_soil();
    ptr->cal_kappa();

    std::vector<double> T_temp = ptr->T;

    // update node temerature
    for (int i = 1; i < this->N - 1; i++) {
        ptr->T[i] = T_temp[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (T_temp[i + 1] - T_temp[i]) - (ptr->kappa[i] + ptr->kappa[i - 1]) / 2 * (T_temp[i] - T_temp[i - 1])) / pow(this->dx, 2) / ptr->rho_c_eff[i];
    }
}
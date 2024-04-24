#include "FDM_icelens_implicit.hpp"

void FDM_icelens_implicit::solve(Ice_model* ptr) {
    // calculate temperature gradient at each node
    ptr->NablaT();

    // calculate rho_c_eff for each node
    ptr->cal_rho_c_eff3();

    // calculate kappa for each node
    ptr->cal_kappa3();

    // // update element temerature
    // for (int i = 1; i < this->N - 1; i++) {
    //     ptr->T[i] = ptr->T[i] + this->dt * ((ptr->kappa[i + 1] + ptr->kappa[i]) / 2 * (ptr->T[i + 1] - ptr->T[i]) - (ptr->kappa[i - 1] + ptr->kappa[i]) / 2 * (ptr->T[i] - ptr->T[i - 1])) / pow(this->dx, 2) / ptr->rho_c_eff[i];
    // }

    arma::mat K(this->N - 1, this->N - 1);
    arma::vec b(this->N - 1);
    if (ptr->idx_lb.back() > 0) {
        for (int  i = 0; i < this->N - 2; i++) {
            b(i) = ptr->rho_c_eff[i + 1] * ptr->T[i + 1];
        }
        b(0) += ptr->T1 * (ptr->kappa[0] + ptr->kappa[1]) / 2 / pow(this->dx, 2) * this->dt;
        b(this->N - 3) += ptr->T2 * (ptr->kappa[this->N - 1] + ptr->kappa[this->N - 2]) / 2 / pow(this->dx, 2) * this->dt;

        (*(this->K))(0, 0) = ptr->rho_c_eff[1] + (ptr->kappa[0] + ptr->kappa[2] + 2 * ptr->kappa[1]) / 2 / pow(this->dx, 2) * this->dt;
        (*(this->K))(0, 1) = - (ptr->kappa[1] + ptr->kappa[2]) / 2 / pow(this->dx, 2) * this->dt;
        for (int i = 1; i < this->N - 3; i++) {
            (*(this->K))(i, i - 1) = - (ptr->kappa[i] + ptr->kappa[i + 1]) / 2 / pow(this->dx, 2) * this->dt;
            (*(this->K))(i, i) = ptr->rho_c_eff[i + 1] + (ptr->kappa[i] + ptr->kappa[i + 2] + 2 * ptr->kappa[i + 1]) / 2 / pow(this->dx, 2) * this->dt;
            (*(this->K))(i, i + 1) = - (ptr->kappa[i + 1] + ptr->kappa[i + 2]) / 2 / pow(this->dx, 2) * this->dt;
        }
        (*(this->K))(this->N - 3, this->N - 4) = - (ptr->kappa[this->N - 2] + ptr->kappa[this->N - 3]) / 2 / pow(this->dx, 2) * this->dt;
        (*(this->K))(this->N - 3, this->N - 3) = ptr->rho_c_eff[this->N - 2] + (ptr->kappa[this->N - 1] + ptr->kappa[this->N - 3] + 2 * ptr->kappa[this->N - 2]) / 2 / pow(this->dx, 2) * this->dt;
        
        std::vector<double> T_new(this->N - 2);
        for (int  i = 0; i < this->N - 2; i++) {
            T_new[i] = ptr->T[i + 1];
        }
    }

    arma::solve(T_new, this->K, b);
    
    for (int  i = 1; i < this->N - 1; i++) {
        ptr->T[i] = T_new[i - 1];
    }   

    // locate Tf
    for (int i = 0; i < this->N - 1; i++) {
        if (((ptr->Tf - 273.15) > ptr->T[i]) && ((ptr->Tf - 273.15) < ptr->T[i + 1])) {
            ptr->xf = ptr->x[i] + (ptr->x[i + 1] - ptr->x[i]) * ((ptr->Tf - 273.15) - ptr->T[i]) / (ptr->T[i + 1] - ptr->T[i]);
            break;
        }
    }

    // get idx during iterations
    ptr->idx = ptr->find_index(ptr->xf);   

    // update time
    ptr->time += this->dt;

    return;
}
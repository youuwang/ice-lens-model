// class Ice_model is a model class including all three phases
// It is a derived class from the abstract class FDM
#ifndef Ice_modelHEADERDEF
#define Ice_modelHEADERDEF

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <algorithm>

using namespace std;

class FDM_ff;
class FDM_solver;
class FDM_icelens;
class Ice_model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
    Ice_model() = default;
    Ice_model(std::string filename) {
        read_input(filename);
        this->Tf = this->Tm - 2 * this->gamma_sl * this->Tm / this->rho_i / this->L / this->Rp;
        this->dx = this->d / (this->N - 1);
        this->dt = pow(dx, 2);
        this->rho_c_eff.assign(this->N, 0);
        this->nablaT.assign(this->N, 0);
        this->kappa_s.assign(this->N, 0);
        this->kappa.assign(this->N, 0);
        this->label.assign(this->N, 1);
        for (int i = 0; i < this->N; i++) {
          this->x.push_back(i * this->dx);
          this->phi_n.push_back(this->phi);
        }
        for (int i = 0; i < this->N; i++) {
          this->T.push_back(this->T01 + (this->T02 - this->T01) / this->d * this->x[i]);
        }
        this->time = 0;
        this->flag = 0;
        this->N_cyc = 1;
        this->time_output = 0;
        this->time_temp = 0;
        this->x_dry = this->d;
        this->water_supply = 1;
    };
    ~Ice_model() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! Solve temperature profile and output results
    void solve(ofstream& output_T, ofstream& output_ice_lens, ofstream& label_time, ofstream& crack_state_time);

    //!read input from .inp file
    void read_input(std::string filename); //read input from .inp file

    //! find index of the element for xx
    int find_index(double xx); 

    //! Calculate temperature gradient at each node
    void NablaT();

  /* ------------------------------------------------------------------------ */
  /* Friend classes                                                           */
  /* ------------------------------------------------------------------------ */ 
public:
    //! delaration of FDM_solver as a friend class to enable calling of protected members in the solver class
    friend class FDM_solver;

    //! delaration of FDM_ff as a friend class to enable calling of protected members in the solver class
    friend class FDM_ff; 

    //! delaration of FDM_icelens as a friend class to enable calling of protected members in the solver class
    friend class FDM_icelens; 


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
    //! soil particle density (kg/m3)
    double rho_s;
    //! Specific heat capacity of soil particles (J/KgK)
    double cp_s;
    //! cohesion of soils (Pa) 
    double c; 
    //! permeability of saturated soils
    double k0;
    //! Soil porosity
    double phi;
    //! Characteristic pore throat size (m)
    double Rp;
    //! Length of the rod (m)
    double d;
    //! Number of grid points
    int N;
    //! Initial temperature on the left (degree Celsius)
    double T01;
    //! Initial temperature on the right (degree Celsius)
    double T02;
    //! Boundary temperature (degree Celsius)
    double T1;
    //! Boundary temperature (degree Celsius)
    double T2;
    //! Rate of cooling or heating at the boundary
    double v_ch;
    //! number of cycles for boundary condition
    int N_cyc;
    //! total number of cycles
    int cyc_tot;
    //! grid size 
    double dx;
    //! time step size
    double dt;
    //! Store coordinates for grid points
    std::vector<double> x;
    //! Store the node temperature
    std::vector<double> T;
    //! store rho_c_eff at each node
    std::vector<double> rho_c_eff;
    //! Thermal conductivity (w/mK)
    std::vector<double> kappa;
    //! Thermal conductivity of soil particles (w/mK)
    std::vector<double> kappa_s;
    //! current position of 0 degree Celsius
    std::vector<double> x0;
    //! speed of 0 degree Celsius
    std::vector<double> v0;
    //! Temperature gradient at grid points
    std::vector<double> nablaT;
    //! lower boundary of ice lens
    std::vector<double> x_lb;
    //! upper boundary of ice lens 
    std::vector<double> x_ub;
    //! index of upper boundary of ice lens
    std::vector<double> idx_ub;
    //! index of lower boundary of ice lens
    std::vector<double> idx_lb;
    //! growth rate of current active ice lens
    double vl;
    //! temperature at lower boundary of current active ice lens
    std::vector<double> Tl;

    //! temperature at upper boundary of ice lenses
    std::vector<double> T_ub;
    
    //! number of 0 degree Celsius
    int N_0;

    //! Boundary for frozen fringe
    double xf;
    //! active x0
    double x0_active;
    //! Temperature at frozen fringe boundary
    double Tf;
    //! drying front
    double x_dry;
    //! time counter
    double time;
    //! indicator for end of cooling or heating
    double time_temp;
    //! indicator for output
    double time_output;
    //! vector of time
    std::vector<double> time_t1;
    std::vector<double> time_t2;
    std::vector<double> time_t3;  
    //! vector of x_ub
    std::vector<double> x_ub_t;
    //! vector of x_lb
    std::vector<double> x_lb_t;
    //! number of ice lenses at different time
    std::vector<int> N_lens;
      //! vector for T_ub
    std::vector<double> T_ub_t;
    //! vector for Tl
    std::vector<double> Tl_t;
    //! vector for xf in phase 2  
    std::vector<double> xf_t2;
    //! vector for xf in phase 3
    std::vector<double> xf_t3;
    //! vector for x_dry
    std::vector<double> x_dry_t;
    //! vector for active ice lens
    std::vector<int> ice_active_t;
    //! Index for the element where Tf is located
    int idx;
    //! Position of the first ice lens
    double lens_ini;

    //! label for each node (1-frozen fringe or saturated medium; 2-pure ice; 3-pure water)
    std::vector<int> label;

    //! porosity at each node (1-phi; 2,3-1) 
    std::vector<double> phi_n;

    //!indicator for heating or cooling process 0-cooling 1-heating
    int flag;
    //! indicator for water supply 0-no water supply 1- have water supply
    int water_supply;

    //! index of element where 0 locates
    std::vector<int> idx_0;

    //! index of active ice lens
    int ice_active;

    //! crack state (0 - pure ice; 1 - half ice half water; 2 - pure water)
    std::vector<int> crack_state;

    //! bool checking whether solve3() is active
    bool solve3_active;

    //! bool checking whether solve2() is active
    bool solve2_active;

    //! Solver pointer of the base class which can also points to any derived class 
    std::shared_ptr<FDM_solver> solver_ptr;

    //! Surface energy at ice-water interface (J/m2)
    const double gamma_sl = 0.03;
    //! Latent heat of fusion for ice (J/Kg)
    const double L = 3.3e5;
    //! Ice density
    const double rho_i = 0.92e3;
    //! Water density
    const double rho_l = 1e3;
    //! Bulk melting temperature for water (K)
    const double Tm = 273.15;
    //! Dynamic viscosity of water(kg/(m*h))
    const double mu = 5e-7;
    //! Specific heat capacity of water (J/KgK)
    const double cp_l = 4.182e3;
    //! Specific heat capacity of ice (J/KgK)
    const double cp_i = 2.108e3;
    //! Thermal conductivity of ice (w/mK)
    const double kappa_i = 7.92e3;
    //! Thermal conductivity of liquid water
    const double kappa_l = 2.052e3;

protected:
    //! solve function for phase 1 (only liquid and soil particles)
    void solve1(ofstream& output_T);
    //! solve function for phase 2 (saturated soil and frozen fringe)
    void solve2(ofstream& output_T);
    //! solve function for phase 3 (ice lens involved following phase 2)
    void solve3(ofstream& output_T, ofstream& label_time, ofstream& crack_state_time);
    //! calculate kappa of soil under different temperatures
    void kappa_soil();
    //! function to calculate kappa of soil particles
    double kappa_soil_fun(double x);
    //! calculate derivative of kappa_soil_fun()
    double kappa_fun_der(double x);
    //! calculate kappa gradient in frozen fringe
    double kappa_gradient(double x);
    //! calculate rho_c_eff for each element
    void cal_rho_c_eff();
    //! calculate kappa for each node
    void cal_kappa();
    //! integral in thermomolecular force
    double FT_int(double x);
    //! integrated function in FT_int
    double FT_int_fun(double x);
    //! integral in expression of hydrodynamic force
    double F_mu_int(double x);
    //! integrated function in F_mu_int
    double F_mu_int_fun(double x);
    //! prescribed ice saturation model (no ice lens)
    std::vector<double> ice_sat(); 
    //! calculate ice saturation for a certain temperature
    double ice_sat_fun(double x);
    //! take derivative of ice saturation function
    double ice_sat_der(double x);
    //! locate ice lens above x0
    int ice_lens_active();


};
#endif
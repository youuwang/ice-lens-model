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
#include <numeric>

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
    Ice_model(std::string filename, std::string odir) {
        read_input(filename);
        this->out_dir = odir;
        this->Tf = this->Tm - 2 * this->gamma_sl * this->Tm / this->rho_i / this->L / this->Rp;
        this->dx.assign(this->N, this->d / (this->N - 1));
        this->dx[0] = 0;
        this->dt = pow(this->d / (this->N - 1), 2);
        this->rho_c_eff.assign(this->N, 0);
        this->nablaT.assign(this->N, 0);
        this->kappa.assign(this->N, 0);
        this->label.assign(this->N, 1);
        this->x_0.assign(this->N, 0);

        std::partial_sum(this->dx.begin(), this->dx.end(), this->x_0.begin());
        this->x = this->x_0;

        for (int i = 0; i < this->N; i++) {
            this->T.push_back(this->T01 + (this->T02 - this->T01) / this->d * this->x[i]);
        }
        this->v = 0;
        this->time = 0;
        this->time_output = 0;
        this->x_dry = this->d;
        this->mu = 5e-7 / this->time_sf;
        this->kappa_i = 7.92e3 * this->time_sf;
        this->kappa_l = 2.052e3 * this->time_sf;
        this->kappa_air = 0.05 * this->time_sf;
        this->kappa_s = 2.16e4 * this->time_sf;
        this->kappa_s_bulk = pow(this->kappa_air, this->phi) * pow(this->kappa_s, 1 - this->phi);
    };
    ~Ice_model() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! Solve temperature profile and output results
    void solve(ofstream& output_T, ofstream& output_ice_lens, ofstream& x_time);

    //!read input from .inp file
    void read_input(std::string filename); //read input from .inp file

    //! find index of the element for xx
    int find_index(double xx); 

    //! Calculate temperature gradient at each node
    void NablaT();

    //! Gaussian quadrature
    double gauss(double (Ice_model::*f)(double x), double a);

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

    //! output directory
    std::string out_dir;

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
    //! Grain size (m)
    double R;
    //! Length of the rod (m)
    double d;
    //! Number of grid points
    int N;
    //! Initial temperature (degree Celsius) on the left
    double T01;
    //! Initial temperature (degree Celsius) on the right
    double T02;
    //! Boundary temperature (degree Celsius)
    double T1;
    //! Boundary temperature (degree Celsius)
    double T2;
    //! periodicity of cooling 
    double v_cd;
    //! grid size 
    std::vector<double> dx;
    //! time step size
    double dt;
    //! Store coordinates for grid points
    std::vector<double> x;
    //! Store initial coordinates for grid points
    std::vector<double> x_0;
    //! Store the node temperature
    std::vector<double> T;
    //! store rho_c_eff at each node
    std::vector<double> rho_c_eff;
    //! Thermal conductivity (w/mK)
    std::vector<double> kappa;
    //! Thermal conductivity of soil particles (w/mK)
    double kappa_s;
    //! Kappa of bulk sandstone
    double kappa_s_bulk;
    //! Current position of 0 degree Celsius
    double x0;
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
    //! temperature at upper boundary of current active ice lens
    std::vector<double> T_ub;

    //! Freezing rate
    double v;
    //! Boundary for frozen fringe
    double xf;
    //! Temperature at frozen fringe boundary
    double Tf;
    //! drying front
    double x_dry;
    //! time counter
    double time;
    //! indicator for output
    double time_output;
    //! Dynamic viscosity of water (Pa h)
    double mu;
    //! Thermal conductivity of ice (w/mK)
    double kappa_i;
    //! Thermal conductivity of liquid water
    double kappa_l;
    //! Thermal conductivity of air
    double kappa_air;

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
    //! vector for x in phase 3
    std::vector<double> x_t3;
    //! vector for x_dry
    std::vector<double> x_dry_t;
    //! Index for the element where frozen fringe boundary is located
    int idx;
    //! Position of the first ice lens
    double lens_ini;

    //! labels for nodes 1 - saturated medium or frozen fringe; 2 - cracks; 3 - dry material
    std::vector<int> label;

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
    //! Specific heat capacity of water (J/KgK)
    const double cp_l = 4.182e3;
    //! Specific heat capacity of ice (J/KgK)
    const double cp_i = 2.108e3;
    //! time scaling factor
    const double time_sf = 1;

protected:
    //! solve function for phase 1 (only liquid and soil particles)
    void solve1(ofstream& output_T);
    //! solve function for phase 2 (saturated soil and frozen fringe)
    void solve2(ofstream& output_T);
    //! solve function for phase 3 (ice lens involved following phase 2)
    void solve3(ofstream& output_T, ofstream& x_time);
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
    //! calculate vl
    double vl_fun();
    //! update grid
    void update_grid();


};
#endif
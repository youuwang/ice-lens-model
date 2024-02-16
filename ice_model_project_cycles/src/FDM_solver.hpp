// FDM_without_s is the class for solving heat conduction equation using FDM without phase change considered

#ifndef FDM_solverHEADERDEF
#define FDM_solverHEADERDEF

#include "Ice_model.hpp"
// #include "matrix.hpp"
// #include "newton_raphson_solver.hpp"
#include <memory>
#include <cassert>

using namespace std;
class Ice_model;
class FDM_solver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
    FDM_solver() = default;
    FDM_solver(double d, int N) {
        this->d = d;
        this->N = N;
        this->dx = d / (N - 1);
        this->dt = pow(this->dx, 2);
    };
    ~FDM_solver() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! solve heat conduction equation when there is no phase change
    virtual void solve(Ice_model* ptr);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
    //! number of nodes
    int N;
    //! length of the domain
    double d;
    //! grid size
    double dx;
    //! time step size
    double dt;
};
#endif

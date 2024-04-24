// FDM_icelens_without_s is the class for solving heat conduction equation using FDM after ice lens is initiated
// while considers no heat source induced by phase change (simplification)

#ifndef FDM_icelens_implicitHEADERDEF
#define FDM_icelens_implicitHEADERDEF

#include "FDM_solver.hpp"
#include <armadillo>

using namespace arma;
class Ice_model;
class FDM_icelens_implicit : public FDM_solver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
    FDM_icelens_implicit() = default;
    FDM_icelens_implicit(double d, int N) : FDM_solver(d, N) {};
    ~FDM_icelens_implicit() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! solve heat conduction equation when there is ice lens formed
    void solve(Ice_model* ptr);

};
#endif

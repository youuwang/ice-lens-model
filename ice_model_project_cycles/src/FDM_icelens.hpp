// FDM_icelens_without_s is the class for solving heat conduction equation using FDM after ice lens is initiated
// while considers no heat source induced by phase change (simplification)

#ifndef FDM_icelensHEADERDEF
#define FDM_icelensHEADERDEF

#include "FDM_solver.hpp"

using namespace std;
class Ice_model;
class FDM_icelens : public FDM_solver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
    FDM_icelens() = default;
    FDM_icelens(double d, int N) : FDM_solver(d, N) {};
    ~FDM_icelens() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! solve heat conduction equation when there is ice lens formed
    void solve(Ice_model* ptr);

};
#endif

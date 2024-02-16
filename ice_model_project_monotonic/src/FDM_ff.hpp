// FDM_with_s is the class for solving heat conduction equation using FDM with phase change considered
// namely heat source term is added due to freezing

#ifndef FDM_ffHEADERDEF
#define FDM_ffHEADERDEF

#include "FDM_solver.hpp"

using namespace std;

class Ice_model;
class FDM_ff : public FDM_solver{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
    FDM_ff() = default;
    FDM_ff(double d, int N) : FDM_solver(d, N)  {};
    ~FDM_ff() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
public:
    //! solve heat conduction equation when there is frozen fringe but no ice lens formed
    void solve(Ice_model* ptr);
};
#endif

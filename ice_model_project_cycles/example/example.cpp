#include <iostream>
#include <fstream>
#include "Ice_model.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    //! declare input file
    std::string filename;
    std::string odir;

    if (argc < 3) {
        std::cerr << "Not enough arguments:" << " ./example <sim-name> <output-directory>" << std::endl;
        return 0;
    }
    else {
        filename = argv[1];
        odir = argv[2];
    }

    //! Create objects for output data files (temperature profile and ice lens information) in directory for results
    ofstream output_T(odir + "/" + "T_profile.dat");
    ofstream output_ice_lens(odir + "/" + "Ice_lens.dat");
    ofstream output_crack_state(odir + "/" + "crack_state_time.dat");
    ofstream label_time(odir + "/" + "label_time.dat");

    //! Create an object of Ice_model class and call solve()
    Ice_model wall_v(filename);
    wall_v.solve(output_T, output_ice_lens, label_time, output_crack_state);
    
    output_T.close();
    output_ice_lens.close();
    label_time.close();
    output_crack_state.close();

    return 0;
}


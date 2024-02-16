#include <iostream>
#include <chrono>
#include <fstream>
#include "Ice_model.hpp"
#include <cstdio>
using namespace std;

// int main(int argc, char* argv[]) {
//     // Record the start time
//     auto start_time = std::chrono::high_resolution_clock::now();
    
//     //! declare input file
//     std::string filename;
//     std::string odir;

//     if (argc < 3) {
//         std::cerr << "Not enough arguments:" << " ./example <sim-name> <output-directory>" << std::endl;
//         return 0;
//     }
//     else {
//         filename = argv[1];
//         odir = argv[2];
//     }

//     //! Create objects for output data files (temperature profile and ice lens information) in directory for results
//     ofstream output_T(odir + "/" + "T_profile.dat");
//     ofstream output_ice_lens(odir + "/" + "Ice_lens.dat");

//     //! Create an object of Ice_model class and call solve()
//     Ice_model wall_v(filename);
//     wall_v.solve(output_T, output_ice_lens);
    
//     output_T.close();
//     output_ice_lens.close();

//     // Record the end time
//     auto end_time = std::chrono::high_resolution_clock::now();

//     // Calculate the duration in milliseconds
//     auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

//     // Print the duration
//     std::cout << "Running time: " << duration << " seconds" << std::endl;

//     return 0;
// }

int main(int argc, char* argv[]) {
    //! declare input file
    std::string filename;
    std::string odir;

    if (argc < 3) {
        std::cerr << "Not enough arguments:" << " ./example <sim-name> <output-dir>" << std::endl;
        return 0;
    }
    else {
        filename = argv[1];
        odir = argv[2];
    }

    //! Create objects for output data files (temperature profile and ice lens information) in directory for results
    ofstream output_T(odir + "/" + "T_profile.dat");
    ofstream output_ice_lens(odir + "/" + "Ice_lens.dat");

    //! Create an object of Ice_model class and call solve()
    Ice_model wall_v(filename, odir);
    wall_v.solve(output_T, output_ice_lens);
    
    output_T.close();
    output_ice_lens.close();

    return 0;
}
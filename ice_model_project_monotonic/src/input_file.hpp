#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// accept a filename string and a section string and return a container filled with strings here called lines
std::vector<std::string> read_input_section(std::string filename_input, const std::string section_input)
{
    // define structure to store data from the section line by line
    // this should be soemthing like a "list of strings"
    std::vector<std::string> lines;

    // define the section you are looking for: name it "section"
    const std::string section = section_input;

    // define dimensions (for now fixed, should be read as well)
    // int dim = 2;

    // ----------------------------------------------------------
    // inspired by code from last week:
    // read file and store each line from the section of interest
    // ----------------------------------------------------------
    
    // character used to designate comment in input file
    char comment_char = '#';

    // input file name
    std::string ifname = filename_input + ".inp";
    
    // open the file with check
    std::ifstream read_file(ifname);
    assert(read_file.is_open());

    // loop over file as long as there is content
    bool in_section = false;
    while(!read_file.eof()) {

        // read the line
        std::string line;
        std::getline(read_file, line);

        // remove comments
        size_t found_comment = line.find_first_of(comment_char);
        if (found_comment != std::string::npos)
        line = line.substr(0,found_comment);
        if (line.empty())
        continue;

        // remove whitespace
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        // check if start or end of section
        std::stringstream sstr(line);
        std::string keyword;
        sstr >> keyword;
        if (keyword == "$section") {
            // section stop
            if (in_section) {
                in_section = false;
            }
            else {
                sstr >> keyword;
                if (keyword == section) {
                    // it is this section
                    in_section = true;
                    continue;
                }
            }
        }

        // if not in section, ignore this line
        if (!in_section) {
            continue;
        }

        // it is part of this section -> store line
        lines.push_back(line);
    }
    read_file.close();

    return lines;
}

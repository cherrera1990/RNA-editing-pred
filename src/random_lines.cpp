/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes a list of files as input and outputs a randomised set of lines from all the input files. There are options to set the number of random files to output, to set groups of a number of consecutive lines
 * to be considered as a single entry and to consider the first line in each file as a header to be respected in the output file.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cstdlib>
using namespace std;

//the main function of the program
int main (int argc, char* argv[]) {
    //implementing -h option
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << " list_of_input_files [options]" << endl;
            cout << "Options: -n followed by the total number of random lines or groups of lines   -g followed by the numer of lines to be grouped    -h if there's a header (should be the same in all the input files)" << endl;
            return 0;
        }
    }
    //initialising the default values for the options
    int n = 1000;
    int g = 1;
    bool h = false;
    //here we'll store all the input file paths
    vector<string> input_files;
    //checking for the options and input file paths
    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        //the -h option takes the first line in all files and doesn't include it in the random pool. The first line of the last input file will be added as header to the output
        if (arg.compare("-h") == 0) h = true;
        //the -n option defines the number of lines or groups of lines to be selected from the random pool
        else if (arg.compare("-n") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << " list_of_input_files [options]" << endl;
                cerr << "Options: -n followed by the total number of random lines or groups of lines   -g followed by the numer of lines to be grouped    -h if there's a header (should be the same in all the input files)" << endl;
                return 1;
            }
            else n = atoi(argv[i]);
        }
        //the -g option defines the number of consecutive lines to be used as a group when selecting from the random pool
        else if (arg.compare("-g") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << " list_of_input_files [options]" << endl;
                cerr << "Options: -n followed by the total number of random lines or groups of lines   -g followed by the numer of lines to be grouped    -h if there's a header (should be the same in all the input files)" << endl;
                return 1;
            }
            else g = atoi(argv[i]);
        }
        //if an argument is not an option, it is added to the input file paths
        else {
            input_files.push_back(arg);
        }
    }
    //checking for at least one input file path
    if (input_files.size() == 0) {
        cerr << "Usage: " << argv[0] << " list_of_input_files [options]" << endl;
        cerr << "Options: -n followed by the total number of random lines or groups of lines   -g followed by the numer of lines to be grouped    -h if there's a header (should be the same in all the input files)" << endl;
        return 1;
    }
    //initialising variables
    ifstream input_file;
    string header;
    //in this vector we will store the random pool of groups of lines
    vector<vector<string> *> line_groups;
    //we create the first line group, we use pointers to vectors to avoid copying whole vectors when saving them into the main vector
    vector<string> *line_group = new vector<string> ();
    //reading all the input files line by line
    string line;
    for (int i = 0; i < input_files.size(); i++) {
        input_file.open(input_files[i].c_str());
        if (not input_file.is_open()) {
            cerr << "could not open the input file " << input_files[i] << endl;
            return 1;
        }
        //if the -h option is active, reading the header line
        if (h) {
            if (not getline(input_file, header)) {
                cerr << "The input file " << input_files[i] << "is empty and does not contain a header" << endl;
                return 1;
            }
        }
        //reading the lines and saving them into the groups
        while (getline(input_file, line)) {
            line_group->push_back(line);
            //when we reach the indicated group size, we save the group and create a new one
            if (line_group->size() == g) {
                line_groups.push_back(line_group);
                line_group = new vector<string> ();
            }
        }
        //if the last group is incomplete, we ignore it and print the lines
        if (not line_group->empty()) cerr << "warning: the file " << input_files[i] << " doesn't have a number of lines that is multiple of the grouping. Ignoring the last incomplete group of lines" << endl;;
        input_file.close();
    }
    //if the -h option is active, we print the header
    if (h) cout << header << endl;
    //we initialise the counter and the random seed
    int l = 0;
    srand(time(0));
    //here we will store the position in the vector of the groups that have been selected in order to avoid repetitions
    vector<int> selected;
    //we select the groups of lines until we reach the desired number or have selected the whole pool
    while (l < n and selected.size() != line_groups.size()) {
        //randomly selecting the position of the group in the vector
        int r = rand()%line_groups.size();
        //checking if the selected position has been previously sellected
        bool found = false;
        int s = 0;
        while (s < selected.size() and not found) {
            if (selected[s++] == r) found = true;
        }
        //if the position has been previously selected, repeat the selection and check until we select a new position
        while (found) {
            r = rand()%line_groups.size();
            found = false;
            s = 0;
            while (s < selected.size() and not found) {
                if (selected[s++] == r) found = true;
            }
        }
        //printing the selected group of lines
        for (int i = 0; i < line_groups[r]->size(); i++) {
            cout << line_groups[r]->at(i) << endl;
        }
        //adding the selected position to the vector of selected positions
        selected.push_back(r);
        l++;
    }
    //printing a warning if there were not enough groups to get to the required number
    if (l != n) cerr << "Warning: there are fewer line groups than expected" << endl;
}

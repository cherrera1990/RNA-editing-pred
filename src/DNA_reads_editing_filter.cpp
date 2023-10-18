/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input file two files with vcf style variants annotation, the first one should be RNA seq variants and the second, genomic variants. This program outputs those
 * variants from the first file that don't appear in the second file. It is mainlyu intendedd to be used to filter out the polymorphisms and get only the RNA editign file.
*/
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <stdlib.h>
using namespace std;

//split a string into a vector using a character as a delimitator, not original code
template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

//split a string into a vector using a character as a delimitator, not original code
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

//a comparer for strings in function version as required for the map container
struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};


//the main function that will get called when the program is run
int main (int argc, char* argv[]) {
    //implementing -h option
    if (argc == 2) {
    string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "RNA_variants_file genomic_variants_file" << endl;
            return 0;
        }
    }
    //checking for the correct number of arguments
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << "RNA_variants_file genomic_variants_file" << endl;
        return 1;
    }
    // reading the paths of the input files from the two first arguments
    const string RNA_in_path(argv[1]);
    const string DNA_in_path(argv[2]);
    
    ifstream RNA_in;
    ifstream DNA_in;
    //opening the DNA variants file
    DNA_in.open(DNA_in_path.c_str());
    if (not DNA_in.is_open()) {
        cerr << "could not open genomic variants file" << endl;
        return 1;
    }
    //in this map, we will save all the variants for each chromosome or scaffold
    map<string, vector<int>, comparer> DNA_vars;
    //reading the DNA variants file line by line
    string DNA_line;
    while (getline(DNA_in, DNA_line)) {
        //checking if the line corresponds to a variant annotation
        if (DNA_line[0] != '#') {
            //parsing the line to obtain the chromosome or scaffold and the variant position
            vector<string> line_parsed = split(DNA_line, '\t');
            string scaf_id = line_parsed[0];
            int pos = atoi(line_parsed[1].c_str());
            //addomg the position to the corresponding entry
            map<string, vector<int> >::iterator DNA_it = DNA_vars.insert(pair<string, vector<int> > (scaf_id, vector<int> ())).first;
            DNA_it->second.push_back(pos);
        }
    }
    DNA_in.close();
    //opening the RNA variants file
    RNA_in.open(RNA_in_path.c_str());
    if (not RNA_in.is_open()) {
        cerr << "could not open RNA variants file" << endl;
        return 1;
    }
    //reading the RNA variants file line by line
    string line;
    while (getline(RNA_in, line)) {
        if (line[0] != 'c') {
            //parsing the line to obtain the scaffold or chromosome and the variant position
            vector<string> line_parsed = split (line, '\t');
            string scaf_id = line_parsed[0];
            int pos = atoi(line_parsed[1].c_str());
            //looking for the variant in the saved DNA variants and writing the line if it isn't found in the DNA variants
            map<string, vector<int> >::iterator DNA_it = DNA_vars.find(scaf_id);
            if (DNA_it == DNA_vars.end()) cout << line << endl;
            else {
                bool found = false;
                int i = 0;
                while (not found and i < DNA_it->second.size()) {
                        found = (pos == DNA_it->second[i++]);
                }
                if (not found) cout << line << endl;
            }
        }
    }
    DNA_in.close();
}

/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a gff file and outputs the length of the whole locus of a gene.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <map>
using namespace std;

//an auxiliary function used in the main split funcion, not original code
template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

//a function to split a string into a vector of strings using a delimiter character, not original code
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

//a string comparer to use within maps
struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};

//the main function that will get called when the program is run//the main function that will get called when the program is run
int main (int argc, char* argv[]) {
    //implementing -h option
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << " gff_file" << endl;
            return 0;
        }
    }
    //checking for the correct number of arguments
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " gtf_file" << endl;
        return 1;
    }
    // reading the path of the input file from the first argument
    const string gtf_in_path(argv[1]);
    //opening the gff file
    ifstream gtf_in;
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open transcript reference file" << endl;
        return 1;
    }
    //reading the gff file line by line
    string gtf_line;
    while (getline(gtf_in, gtf_line)) {
        //checking if the line corresponds to an annotation
        if (gtf_line[0] != '#') {
            //parsing the line by tabs
            vector<string> gtf_parsed = split(gtf_line, '\t');
            //checking if the line corresponds to a gene annotation
            bool is_gene = (gtf_parsed[2].compare("gene") == 0);
            if (is_gene) {
                //getting the scaffold id or chromosome
                string scaf_id = gtf_parsed[0];
                //parsing the subfields
                string gene_id;
                vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ';');
                vector<string> gtf_subfield_parsed;
                //checking for the field with the gene id
                for (int i = 0; i < gtf_last_field_parsed.size(); i++) {
                    gtf_subfield_parsed = split(gtf_last_field_parsed[i], '=');
                    if (gtf_subfield_parsed[0].compare("gene_id") == 0 or gtf_subfield_parsed[0].compare("ID") == 0) gene_id = gtf_subfield_parsed[1];
                }
                //using the original gene id as gene name
                string gene_name = gene_id;
                //adding the scaf_id to the gene_id
                gene_id += "_" + scaf_id;
                //getting the annotated start and end positions
                int start, end;
                start = atoi(gtf_parsed[3].c_str());
                end = atoi(gtf_parsed[4].c_str());
                //calculating the length
                int length = end - start + 1;
                //printing the gene name, gene id and the length
                cout << gene_name << "\t" << gene_id << "\t" << length << endl;
            }
        }
    }
    gtf_in.close();
}

/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a gtf file and outputs the length of the whole locus of a gene, taking into consideration the minimum start and maximum end positions of all the annotated transcripts of each gene.
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

//struct containing the relevant information about the gene annotation
typedef struct {
    string gene_id;
    string gene_name;
    int gtf_start;
    int gtf_end;
    bool positive_strand;
    vector<string> editing_lines;
    string seq;
} gtf_info;


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
            cout << "Usage: " << argv[0] << " gtf_file" << endl;
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
    //opening the gtf file
    ifstream gtf_in;
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open transcript reference file" << endl;
        return 1;
    }
    //in this map we will store the output channels for each gene
    map<string, vector<gtf_info>, comparer> amphi_gtf_info;
    //this map will hold the gene ids and the number of different annotation found for those gene ids
    map<string, int> gene_alts;
    //reading the gtf file line by line
    string gtf_line;
    while (getline(gtf_in, gtf_line)) {
        //checking if the line corresponds to an annotation
        if (gtf_line[0] != '#') {
            //parsing the line by tabs
            vector<string> gtf_parsed = split(gtf_line, '\t');
            //checking if the line corresponds to a transcript annotation
            bool is_gene = (gtf_parsed[2].compare("transcript") == 0);
            if (is_gene) {
                //getting the scaffold id or chromosome
                string scaf_id = gtf_parsed[0];
                //creating the gtf info entry
                gtf_info n_info;
                //parsing the subfields
                vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ';');
                vector<string> gtf_subfield_parsed;
                //searching for the gene_id in the subfields
                for (int i = 0; i < gtf_last_field_parsed.size(); i++) {
                    if (gtf_last_field_parsed[i][0] == ' ') gtf_last_field_parsed[i].erase(0,1);
                    gtf_subfield_parsed = split(gtf_last_field_parsed[i], ' ');
                    if (gtf_subfield_parsed[0].compare("gene_id") == 0) n_info.gene_id = gtf_subfield_parsed[1];
                }
                //removing the quotation marks
                n_info.gene_id.erase(n_info.gene_id.length() - 1, 1);
                n_info.gene_id.erase(0, 1);
                //the gene_name is the original gene_id
                n_info.gene_name = n_info.gene_id;
                //for the gene_id we use the gene_id concatenated to the scaf_id
                n_info.gene_id += "_" + scaf_id;
                //getting the annotated start and end positions and the strand
                n_info.gtf_start = atoi(gtf_parsed[3].c_str());
                n_info.gtf_end = atoi(gtf_parsed[4].c_str());
                n_info.positive_strand = (gtf_parsed[6][0] == '+');
                //if the current scaffold hasn't been seen, we insert a new entry, otherwise, we get the corresponding entry
                map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.insert(pair<string, vector<gtf_info> > (scaf_id, vector<gtf_info> ())).first;
                //finging if we already have information on the current gene_id
                bool found = false;
                int i = 0;
                while (not found and i < gtf_it->second.size()) {
                    if (n_info.gene_name.compare(gtf_it->second[i].gene_name) == 0) found = true;
                    else i++;
                }
                //if the gene_id already exists, we check if the start or end need to be updated
                if (found) {
                    if (n_info.gtf_start < gtf_it->second[i].gtf_start) gtf_it->second[i].gtf_start = n_info.gtf_start;
                    if (n_info.gtf_end < gtf_it->second[i].gtf_end) gtf_it->second[i].gtf_end = n_info.gtf_end;
                }
                //if the gene_id doesn't exist, we add the corresponding entry
                else {
                    gtf_it->second.push_back(n_info);
                }
            }
        }
    }
    gtf_in.close();
    //we go over all the entries, printing the gene_name, gene_id and calculated premrna length
    map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.begin();
    while (gtf_it != amphi_gtf_info.end()) {
        for (int i = 0; i < gtf_it->second.size(); i++) {
            int start = gtf_it->second[i].gtf_start;
            int end = gtf_it->second[i].gtf_end;
            int length = end - start + 1;
            cout << gtf_it->second[i].gene_name << "\t" << gtf_it->second[i].gene_id << "\t" << length << endl;
        }
        gtf_it++;
    }
}

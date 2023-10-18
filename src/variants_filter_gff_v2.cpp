/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes the variants specified in vcf format (as given by the SAMtools library) and filters the variants, only printing those that correspond
 * to changes from A to G or T to C, which are candidates for A-to-I RNA-editing events. It uses a gff file as an input to find the genes and transcripts that 
 * contain the events. Several options can be selected for further filtering criteria, including depth, score and whether the variant is within a coding region
 * or not.
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

//struct to save the relevant info from a the gene annotation
typedef struct {
    string gene_id;
    string trans_id;
    int gtf_start;
    int gtf_end;
    bool positive_strand;
    bool CDS;
} gtf_info;

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
            cout << "Usage: " << argv[0] << "vcf_file gtf_file [options]" << endl;
            cout << "Usage: " << " options: -s score filter      -d depth filter       -v (vcf only output)       -i (ignore variants not presence in gtf file)       -c (ignore variants not presence marked as CDS in gtf file)" << endl;
            return 0;
        }
    }
    // checking if the minimum number of arguments is correct
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << "vcf_file gtf_file [options]" << endl;
        return 1;
        cerr << "Usage: " << " options: -s score filter      -d depth filter       -v (vcf only output)       -i (ignore variants not presence in gtf file)       -c (ignore variants not presence marked as CDS in gtf file)" << endl;
    }
    // reading the paths of the input files from the two first arguments
    const string vcf_in_path(argv[1]);
    const string gtf_in_path(argv[2]);
    //initialising parameters for implementing options
    bool depth_filter = false;
    int min_depth = -1;
    bool score_filter = false;
    int min_score = -1;
    bool vcf_only = false;
    bool ignore_nongene = false;
    bool ignore_noncoding = false;
    // processing options
    for (int i = 3; i < argc; i++) {
        if (string(argv[i]) == "-s") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << "vcf_file gtf_file [options]" << endl;
                cerr << "Usage: " << " options: -s score filter      -d depth filter       -v (vcf only output)       -i (ignore variants not presence in gtf file)       -c (ignore variants not presence marked as CDS in gtf file)" << endl;
                return 1;
            }
            min_score = atoi(argv[i]);
            score_filter = true;
        }
        else if (string(argv[i]) == "-d") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << "vcf_file gtf_file [options]" << endl;
                cerr << "Usage: " << " options: -s score filter      -d depth filter       -v (vcf only output)       -i (ignore variants not presence in gtf file)       -c (ignore variants not presence marked as CDS in gtf file)" << endl;
                return 1;
            }
            min_depth = atoi(argv[i]);
            depth_filter = true;
        }
        else if (string(argv[i]) == "-v") {
            vcf_only = true;
        }
        else if (string(argv[i]) == "-i") {
            ignore_nongene = true;
        }
        else if (string(argv[i]) == "-c") {
            ignore_noncoding = true;
        }
        else {
            cerr << "Usage: " << argv[0] << "vcf_file gtf_file [options]" << endl;
            cerr << "Usage: " << " options: -s score filter      -d depth filter       -v (vcf only output)       -i (ignore variants not presence in gtf file)       -c (ignore variants not presence marked as CDS in gtf file)" << endl;
            return 1;
        }
    }
    //declaring input file streams
    ifstream vcf_in;
    ifstream gtf_in;
    //opening gff input file
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open transcript reference file" << endl;
        return 1;
    }
    //this map will save all the info of each gff line with the scaffold id as the key and a vector of the gtf_info struct
    map<string, vector<gtf_info>, comparer> amphi_gtf_info;
    //reading the gtf file line by line
    string gtf_line;
    while (getline(gtf_in, gtf_line)) {
        //ignoring the lines that don't include an annotation
        if (gtf_line[0] != '#') {
            //splitting the line by tab
            vector<string> gtf_parsed = split(gtf_line, '\t');
            //checking if the line corresponds to a gene annotation
            bool is_gene = (gtf_parsed[2].compare("gene") == 0);
            if (is_gene) {
                //saving all the relevant info from the gtf line into the gtf_info struct
                string scaf_id = gtf_parsed[0];
                gtf_info n_info;
                vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ';');
                vector<string> gff_subfield_parsed;
                for (int i = 0; i < gtf_last_field_parsed.size(); i++) {
                    gff_subfield_parsed = split(gtf_last_field_parsed[i], '=');
                    if (gff_subfield_parsed[0].compare("gene_id") == 0) n_info.gene_id = gff_subfield_parsed[1];
                }
                n_info.trans_id = n_info.gene_id;
                n_info.gtf_start = atoi(gtf_parsed[3].c_str());
                n_info.gtf_end = atoi(gtf_parsed[4].c_str());
                n_info.positive_strand = (gtf_parsed[6][0] == '+');
                n_info.CDS = false;
                //if the map doesn't have the scaffold entry yet, we insert a new entry for the scaffold, if it already exists we retrieve it. 
                map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.insert(pair<string, vector<gtf_info> > (scaf_id, vector<gtf_ino> ())).first;
                // we add the new gtf line in the vector for the scaffold map entry
                gtf_it->second.push_back(n_info);
            }
        }

    }
    //opening the vcf file
    vcf_in.open(vcf_in_path.c_str());
    if (not vcf_in.is_open()) {
        cerr << "could not open vcf file" << endl;
        return 1;
    }
    //reading the vcf file line by line
    string line;
    bool eof = false;
    while (getline(vcf_in, line)) {
        //ignoring header lines starting with '#'
        if (line[0] != '#') {
            //spliting the line by tab
            vector<string> parsed = split (line, '\t');
            //getting the position within the scaffold
            int pos = atoi(parsed[1].c_str());
            //parsing the 8th field in the VCF line, containing diverse information separated by ';'
            vector<string> param_parsed = split (parsed[7], ';');
            //ignoring indels
            if (parsed[7][0] != 'I') {
                //checking the depth and score if the corresponding options were selected
                vector<string> depth_parsed = split (param_parsed[0], '=');
                int depth = atoi(depth_parsed[1].c_str());
                double score = atof(parsed[5].c_str());
                if ((not score_filter or score >= (double) min_score) and (not depth_filter or depth >= min_depth)) {
                    //declaring vectors that will hold the information about the transcripts from gtf that contain the position of the variant in the current vcf line
                    vector<vector<string> > trans_ids;
                    vector<vector<string> > coding_trans_ids;
                    vector<bool> positive_strands;
                    vector<string> gene_ids;
                    map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.find(parsed[0]);
                    //searching for the scaffold in the gtf info map
                    if (gtf_it == amphi_gtf_info.end()) {
                        cerr << "Something went wrong, Unable to find scaffold " << parsed[0] << " in gtf file" << endl; 
                    }
                    else {
                        //we check all the gff info structs that match the scaffold from the current vcf line
                        for (int i = 0; i < gtf_it->second.size(); i++) {
                            //we check if the curent variant is within the positions specified by the gtf file
                            if (pos >= gtf_it->second[i].gtf_start and pos <= gtf_it->second[i].gtf_end) {
                                //if the current variant has only one alternative nucleotide, we check if the change is either A to G or T to C
                                if (parsed[4].length() == 1) {
                                    if (gtf_it->second[i].positive_strand and parsed[3][0] == 'A' and parsed[4][0] == 'G'
                                      or not gtf_it->second[i].positive_strand and parsed[3][0] == 'T' and parsed[4][0] == 'C') {
                                        //we check if we already have found the current gene in the gff info associated to the current variant
                                        int j = 0;
                                        bool gene_found = false;
                                        while (not gene_found and j < gene_ids.size()) {
                                            if (gene_ids[j].compare(gtf_it->second[i].gene_id) == 0) gene_found = true;
                                            else j++;
                                        }
                                        //if we haven't found the current gene for this variant, we add the information about it
                                        if (not gene_found) {
                                            gene_ids.push_back(gtf_it->second[i].gene_id);
                                            positive_strands.push_back(gtf_it->second[i].positive_strand);
                                            trans_ids.push_back(vector<string> ());
                                            coding_trans_ids.push_back(vector<string> ());
                                        }
                                        //we add the transcript corresponding to the current gtf info field associated to the current variant, differentiating between coding and noncoding regions
                                        if (not gtf_it->second[i].CDS) trans_ids[j].push_back(gtf_it->second[i].trans_id);
                                        else coding_trans_ids[j].push_back(gtf_it->second[i].trans_id);
                                    }
                                }
                                //if the current variant has two alternative nucleotides we check if one of the changes corresponds to possible editing
                                else if (parsed[4].length() == 3) {
                                    if (gtf_it->second[i].positive_strand and parsed[3][0] == 'A' and (parsed[4][0] == 'G' or parsed[4][2] == 'G')
                                      or not gtf_it->second[i].positive_strand and parsed[3][0] == 'T' and (parsed[4][0] == 'C' or parsed[4][2] == 'C')) {
                                        //we check if we already have found the current gene in the gff info associated to the current variant
                                        int j = 0;
                                        bool gene_found = false;
                                        while (not gene_found and j < gene_ids.size()) {
                                            if (gene_ids[j].compare(gtf_it->second[i].gene_id) == 0) gene_found = true;
                                            else j++;
                                        }
                                        //if we haven't found the current gene for this variant, we add the information about it
                                        if (not gene_found) {
                                            gene_ids.push_back(gtf_it->second[i].gene_id);
                                            positive_strands.push_back(gtf_it->second[i].positive_strand);
                                            trans_ids.push_back(vector<string> ());
                                            coding_trans_ids.push_back(vector<string> ());
                                        }
                                        //we add the transcript corresponding to the current gtf info field associated to the current variant, differentiating between coding and noncoding regions
                                        if (not gtf_it->second[i].CDS) trans_ids[j].push_back(gtf_it->second[i].trans_id);
                                        else coding_trans_ids[j].push_back(gtf_it->second[i].trans_id);
                                    }
                                }
                            }
                        }
                        //we check if the vcf format output option was selected (we print the original lines, without adding the information about genes and transcripts
                        if (vcf_only) {
                            //we check for different conditions depending on the options selected and print the vcf line if all conditions are met
                            if (not ignore_nongene and gene_ids.size() == 0) {
                                if (parsed[4].length() == 1 and (parsed[3][0] == 'A' and parsed[4][0] == 'G'
                                  or parsed[3][0] == 'T' and parsed[4][0] == 'C')) cout << line << endl;
                                else if (parsed[4].length() == 3 and (parsed[3][0] == 'A' and (parsed[4][0] == 'G' or parsed[4][2] == 'G')
                                  or parsed[3][0] == 'T' and (parsed[4][0] == 'C' or parsed[4][2] == 'C'))) { 
                                    cout << line << endl;
                                }
                            }
                            else if (not ignore_noncoding) {
                                bool trans_found = false;
                                int i = 0;
                                while (not trans_found and i < trans_ids.size()) {
                                    trans_found = trans_ids[i++].size() > 0;
                                }
                                if (trans_found) cout << line << endl;
                            }
                            else {
                                bool trans_found = false;
                                int i = 0;
                                while (not trans_found and i < coding_trans_ids.size()) {
                                    trans_found = coding_trans_ids[i++].size() > 0;
                                }
                                if (trans_found) cout << line << endl;
                            }
                        }
                        //if the vcf output option was not selected we print the lines adding the information about the genes, transcripts and strands that were found associated to the variant
                        else {
                            if (gene_ids.size() > 0) {
                                for (int i = 0; i < gene_ids.size(); i++) {
                                    if (not ignore_noncoding and trans_ids[i].size() > 0) {
                                        cout << line;
                                        if (positive_strands[i]) cout << "\t+\t" << gene_ids[i] << '\t';
                                        else cout << "\t-\t" << gene_ids[i] << '\t';
                                        for (int j = 0; j < trans_ids[i].size() - 1; j++) {
                                            cout << trans_ids[i][j] << ":";
                                        }
                                        cout << trans_ids[i][trans_ids[i].size() - 1] << "\texon" << endl;
                                    }
                                }
                            }
                            if (coding_trans_ids.size() > 0) {
                                for (int i = 0; i < gene_ids.size(); i++) {
                                    if (coding_trans_ids[i].size() > 0) {
                                        cout << line;
                                        if (positive_strands[i]) cout << "\t+\t" << gene_ids[i] << '\t';
                                        else cout << "\t-\t" << gene_ids[i] << '\t';
                                        for (int j = 0; j < coding_trans_ids[i].size() - 1; j++) {
                                            cout << coding_trans_ids[i][j] << ":";
                                        }
                                        cout << coding_trans_ids[i][coding_trans_ids[i].size() - 1] << "\tcoding" << endl;
                                    }
                                }
                            }
                            if (not ignore_nongene and gene_ids.size() == 0) {
                                if (parsed[4].length() == 1 and (parsed[3][0] == 'A' and parsed[4][0] == 'G'
                                  or parsed[3][0] == 'T' and parsed[4][0] == 'C')) cout << line << "\t*\tunknown_gene\tunknown_trans\t*" << endl;
                                else if (parsed[4].length() == 3 and (parsed[3][0] == 'A' and (parsed[4][0] == 'G' or parsed[4][2] == 'G')
                                  or parsed[3][0] == 'T' and (parsed[4][0] == 'C' or parsed[4][2] == 'C'))) {
                                    cout << line << "\t*\tunknown_gene\tunknown_trans\t*" << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
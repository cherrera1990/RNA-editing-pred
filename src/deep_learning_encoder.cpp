/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a file with a list of positions within genes as given by the programs full_premrna_sequence.cpp and non_editing_adenosines.cpp, the secondary structure file with a dot-parenthesis
 * notation as given by linearfold and the secondary structure annotationn file as given by  the structure_features.cpp program and gives for each gene, a series of tracks. First, the title of the sequence is
 * given (as in the secondary structure file), then the first and second tracks are the nucleotide and dot-parenthesis notation of the structure, exactly as given by the input file. The third track has
 * the annotated features corresponding to each position, annotated with a single letter: i for inner loops, b for bulge loops, h for hairpin loops, d for double strand and s for single strand. The fourth and
 * final track annotates the positions listed as input with a 2, other adenosines with a 1 and other nucleotides with a 0. There's an empty line between the tracks of each gene.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <stdlib.h>
#include <map>

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

//struct containing the sequence and all the output channels
typedef struct {
    string seq;
    string str;
    string edit;
    string str_type;
} dl_info;

//a comparer for strings in function version as required for the map container
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
            cout << "Usage: " << argv[0] << " RNA_editing_file structures_file structures_annotation_file" << endl;
            return 0;
        }
    }
    //checking for the correct number of arguments
    if (argc != 4) {
            cerr << "Usage: " << argv[0] << " RNA_editing_file structures_file structures_annotation_file" << endl;
            return 1;
    }
    // reading the paths of the input files from the three first arguments
    const string edit_in_path(argv[1]);
    const string str_in_path(argv[2]);
    const string structures_annot_in_path(argv[3]);
    ifstream edit_in;
    ifstream str_in;
    ifstream str_annot_in;
    //opening the structure input file
    str_in.open(str_in_path.c_str());
    if (not str_in.is_open()) {
        cerr << "could not open RNA structures file" << endl;
        return 1;
    }
    //in this map we will store the output channels for each gene
    map<string, dl_info, comparer> deep_learning_info;
    //declaring variables
    string title_line;
    string seq;
    string str;
    //reading the title line in the structures file
    while (getline(str_in, title_line)) {
        //reading the sequence and structure lines for each gene and giving error if there are no lines to read
        if (not getline(str_in, seq)) {
            cerr << "error " << title_line << " has no sequence" << endl;
            return 1;
        }
        if (not getline(str_in, str)) {
            cerr << "error " << title_line << " has no structure" << endl;
            return 1;
        }
        //creating the struct containing the cnannels for the gene
        dl_info dli;
        //copying the sequence and structure into the struct
        dli.seq = seq;
        dli.str = str;
        //inirialising the editing channel with 0s as non-adenosine
        dli.edit = string(seq.length(), '0');
        dli.str_type = string(seq.length(), '0');
        //looking for adenosines in the sequence and changing the positions to 1 as non-edited adenosine
        for (int i = 0; i < dli.edit.length(); i++) {
            if (seq[i] == 'a' or seq[i] == 'A') dli.edit[i] = '1';
        }
        //processing the title line to obtain the gene ID
        title_line.erase(0, 1);
        string title = title_line;
        //saving the struct into the map
        deep_learning_info.insert(pair<string, dl_info> (title, dli));
    }
    str_in.close();
    //opening the editing annotation file
    edit_in.open(edit_in_path.c_str());
    if (not edit_in.is_open()) {
        cerr << "could not open RNA editing file" << endl;
        return 1;
    }
    //reading the annotation file line by line
    string line;
    while (getline(edit_in, line)) {
        //parsing the line to obtain the gene id and the position inside the gene
        vector<string> parsed = split(line, '\t');
        string gene_id = parsed[0] + "\t" + parsed[1];
        int pos_in_gene = atoi(parsed[2].c_str());
        //finding the struct containing the channels for the particular gene
        map<string, dl_info>::iterator dli_it = deep_learning_info.find(gene_id);
        if (dli_it == deep_learning_info.end()) {
            cerr << "Warning: Something went wrong, Unable to find gene " << gene_id << " in structures file" << endl; 
        }
        else {
            //marking the position in the editing channel as an edited adenosine
            dli_it->second.edit[pos_in_gene - 1] = '2';
        }
    }
    edit_in.close();
    //opening the structure features annotation file
    str_annot_in.open(structures_annot_in_path.c_str());
    if (not str_annot_in.is_open()) {
        cerr << "could not open RNA structures annotation file" << endl;
        return 1;
    }
    //reading the structure features annotation file line by line
    string structures_line;
    map<string, dl_info>::iterator current_structure_it;
    while (getline(str_annot_in, structures_line)) {
        //checking if it is not a start of a new gene annotation
        if (structures_line[0] != '>') {
            //parsing the line using spaces as separators
            vector<string> structures_parsed = split(structures_line, ' ');
            //checking if it is a double strand fragment annotation
            if (structures_parsed[0].compare("DS") == 0) {
                //parsing the startt and end positions of both strands of the DS fragment
                vector<string> p1 = split(structures_parsed[1], '=');
                vector<string> p1_parsed = split(p1[1], ',');
                int start1 = atoi(p1_parsed[0].c_str()) - 1;
                int start2 = atoi(p1_parsed[1].c_str()) - 1;
                vector<string> p2 = split(structures_parsed[2], '=');
                vector<string> p2_parsed = split(p2[1], ',');
                int end1 = atoi(p2_parsed[0].c_str()) - 1;
                int end2 = atoi(p2_parsed[1].c_str()) - 1;
                vector<string> s = split(structures_parsed[3], '=');
                //annotating the corresponding positions with 'd' for double strand
                current_structure_it->second.str_type.replace(start1, end1 - start1 + 1, end1 - start1 + 1, 'd');
                current_structure_it->second.str_type.replace(end2, start2 - end2 + 1, start2 - end2 + 1, 'd');
                
            }
            //checking if it is a single strand fragment annotation
            else if (structures_parsed[0].compare("SSMF") == 0 or structures_parsed[0].compare("SSBF") == 0 or structures_parsed[0].compare("SSEF") == 0) {
                //parsing the start and end positions of the single strand fragment
                vector<string> st = split(structures_parsed[1], '=');
                int start1 = atoi(st[1].c_str()) - 1;
                vector<string> en = split(structures_parsed[2], '=');
                int end1 = atoi(en[1].c_str()) - 1;
                //annotating the corresponding positions with 's' for single strand
                current_structure_it->second.str_type.replace(start1, end1 - start1 + 1, end1 - start1 + 1, 's');
            }
            //checking if it is a hairpin loop annotation
            else if (structures_parsed[0].compare("HL") == 0) {
                //parsing the start and end positions of the hairpin loop
                vector<string> st = split(structures_parsed[1], '=');
                int start1 = atoi(st[1].c_str()) - 1;
                vector<string> en = split(structures_parsed[2], '=');
                int end1 = atoi(en[1].c_str()) - 1;
                //annotating the corresponding positions with 'h' for hairpin loop
                current_structure_it->second.str_type.replace(start1, end1 - start1 + 1, end1 - start1 + 1, 'h');
            }
            //the branching point annotation is not used anymore
            /*else if (structures_parsed[0].compare("BP") == 0) {
            }*/
            //checking if it is an inner loop annotation
            else if (structures_parsed[0].compare("IL") == 0 or structures_parsed[0].compare("ILN") == 0) {
                //parsing the start and end positions of both strands of the inner loop
                bool nick = (structures_parsed[0].compare("ILN") == 0);
                vector<string> st1 = split(structures_parsed[1], '=');
                int start1 = atoi(st1[1].c_str()) - 1;
                vector<string> en1 = split(structures_parsed[2], '=');
                int end1 = atoi(en1[1].c_str()) - 1;
                vector<string> st2 = split(structures_parsed[5], '=');
                int start2 = atoi(st2[1].c_str()) - 1;
                vector<string> en2 = split(structures_parsed[6], '=');
                int end2 = atoi(en2[1].c_str()) - 1;
                //annotating the corresponding positions with 'i' for inner loop
                current_structure_it->second.str_type.replace(start1, end1 - start1 + 1, end1 - start1 + 1, 'i');
                current_structure_it->second.str_type.replace(end2, start2 - end2 + 1, start2 - end2 + 1, 'i');
            }
            //checking if it is an bulge loop annotation
            else if (structures_parsed[0].compare("BL") == 0 or structures_parsed[0].compare("BLN") == 0) {
                //parsing the start and end positions of both strands of the bulge loop
                bool nick = (structures_parsed[0].compare("BLN") == 0);
                vector<string> st1 = split(structures_parsed[1], '=');
                int start1 = atoi(st1[1].c_str()) - 1;
                vector<string> en1 = split(structures_parsed[2], '=');
                int end1 = atoi(en1[1].c_str()) - 1;
                vector<string> st2 = split(structures_parsed[5], '=');
                int start2 = atoi(st2[1].c_str()) - 1;;
                vector<string> en2 = split(structures_parsed[6], '=');
                int end2 = atoi(en2[1].c_str()) - 1;
                //annotating the corresponding positions with 'i' for inner loop
                current_structure_it->second.str_type.replace(start1, end1 - start1 + 1, end1 - start1 + 1, 'b');
                current_structure_it->second.str_type.replace(end2, start2 - end2 + 1, start2 - end2 + 1, 'b');
            }
            //ignoring the stats lines
            else if (structures_parsed[0].compare("DSstats") == 0) {
                
            }
            else if (structures_parsed[0].compare("HLstats") == 0) {
                
            }
            else if (structures_parsed[0].compare("ILstats") == 0) {
                
            }
            else if (structures_parsed[0].compare("BLstats") == 0) {
                
            }
            else cerr << "Warning: line \"" << structures_line << "\" not recognised" << endl; 
        }
        //when we find a new gene annotation, we change parse the gene id and change the current map entry to the new gene
        else {
            structures_line.erase(0, 1);
            current_structure_it = deep_learning_info.find(structures_line);
        }
    }
    str_annot_in.close();
    //we go over all the gene entries and print the title and all the corresponding channels
    for (map<string, dl_info>::iterator dli_it = deep_learning_info.begin(); dli_it != deep_learning_info.end(); dli_it++) {
        cout << dli_it->first << endl;
        cout << dli_it->second.seq << endl;
        cout << dli_it->second.str << endl;
        cout << dli_it->second.str_type << endl;
        cout << dli_it->second.edit << endl;
        cout << endl;
    }
}
/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a file with a list of positions within genes as given by the programs full_premrna_sequence.cpp and non_editing_adenosines.cpp, the secondary structure file with a dot-parenthesis
 * notation as given by linearfold and gives the list of the positions of all adenosines that are not annotated in the input positions. In the current version, the output lines are divided among 40 files
 * (in the future, an option to change this number may be implemented).
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

//struct containing the sequence and the structure and editing channels
typedef struct {
    string seq;
    string str;
    string edit;
} dl_info;

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
            cout << "Usage: " << argv[0] << " RNA_editing_file structures_file output_files_path" << endl;
            return 0;
        }
    }
    //checking for the correct number of arguments
    if (argc != 4) {
            cerr << "Usage: " << argv[0] << " RNA_editing_file structures_file output_files_path" << endl;
            return 1;
    }
    // reading the paths of the input files from the three first arguments
    const string edit_in_path(argv[1]);
    const string str_in_path(argv[2]);
    const string output_file_path(argv[3]);
    //we have two input files and the output divided among 40 files
    ifstream edit_in;
    ifstream str_in;
    ofstream out0, out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12, out13, out14, out15, out16, out17, out18, out19, out20, out21, out22, out23, out24, out25, out26, out27, out28, out29, out30, out31, out32, out33, out34, out35, out36, out37, out38, out39;
    //opening the structure input file
    str_in.open(str_in_path.c_str());
    if (not str_in.is_open()) {
        cerr << "could not open RNA structures file" << endl;
        return 1;
    }
    //in this map we will store the output channels for each gene
    map<string, dl_info, comparer> deep_learning_info;
    vector<string> gene_ids;
    //declaring variables to read lines
    string title_line;
    string seq;
    string str;
    //in this map we will store the output channels for each gene
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
        //inirialising the editing channel with 0s as non-adenosine
        for (int i = 0; i < dli.edit.length(); i++) {
            if (seq[i] == 'a' or seq[i] == 'A') dli.edit[i] = '1';
        }
        title_line.erase(0, 1);
        //processing the title line to obtain the gene ID
        string title = title_line;
        //saving the struct into the map
        deep_learning_info.insert(pair<string, dl_info> (title, dli));
    }
    str_in.close();
    //opening the editing input file and all the output files
    edit_in.open(edit_in_path.c_str());
    if (not edit_in.is_open()) {
        cerr << "could not open RNA editing file" << endl;
        return 1;
    }
    out0.open(output_file_path + "_0", ofstream::trunc);
    if (not out0.is_open()) {
        cerr << "could not open " + output_file_path + "_0 output file" << endl;
        return 1;
    }
    out1.open(output_file_path + "_1", ofstream::trunc);
    if (not out1.is_open()) {
        cerr << "could not open " + output_file_path + "_1 output file" << endl;
        return 1;
    }
    out2.open(output_file_path + "_2", ofstream::trunc);
    if (not out2.is_open()) {
        cerr << "could not open " + output_file_path + "_2 output file" << endl;
        return 1;
    }
    out3.open(output_file_path + "_3", ofstream::trunc);
    if (not out3.is_open()) {
        cerr << "could not open " + output_file_path + "_3 output file" << endl;
        return 1;
    }
    out4.open(output_file_path + "_4", ofstream::trunc);
    if (not out4.is_open()) {
        cerr << "could not open " + output_file_path + "_4 output file" << endl;
        return 1;
    }
    out5.open(output_file_path + "_5", ofstream::trunc);
    if (not out5.is_open()) {
        cerr << "could not open " + output_file_path + "_5 output file" << endl;
        return 1;
    }
    out6.open(output_file_path + "_6", ofstream::trunc);
    if (not out6.is_open()) {
        cerr << "could not open " + output_file_path + "_6 output file" << endl;
        return 1;
    }
    out7.open(output_file_path + "_7", ofstream::trunc);
    if (not out7.is_open()) {
        cerr << "could not open " + output_file_path + "_7 output file" << endl;
        return 1;
    }
    out8.open(output_file_path + "_8", ofstream::trunc);
    if (not out8.is_open()) {
        cerr << "could not open " + output_file_path + "_8 output file" << endl;
        return 1;
    }
    out9.open(output_file_path + "_9", ofstream::trunc);
    if (not out9.is_open()) {
        cerr << "could not open " + output_file_path + "_9 output file" << endl;
        return 1;
    }
    out10.open(output_file_path + "_10", ofstream::trunc);
    if (not out10.is_open()) {
        cerr << "could not open " + output_file_path + "_10 output file" << endl;
        return 1;
    }
    out11.open(output_file_path + "_11", ofstream::trunc);
    if (not out11.is_open()) {
        cerr << "could not open " + output_file_path + "_11 output file" << endl;
        return 1;
    }
    out12.open(output_file_path + "_12", ofstream::trunc);
    if (not out12.is_open()) {
        cerr << "could not open " + output_file_path + "_12 output file" << endl;
        return 1;
    }
    out13.open(output_file_path + "_13", ofstream::trunc);
    if (not out13.is_open()) {
        cerr << "could not open " + output_file_path + "_13 output file" << endl;
        return 1;
    }
    out14.open(output_file_path + "_14", ofstream::trunc);
    if (not out14.is_open()) {
        cerr << "could not open " + output_file_path + "_14 output file" << endl;
        return 1;
    }
    out15.open(output_file_path + "_15", ofstream::trunc);
    if (not out15.is_open()) {
        cerr << "could not open " + output_file_path + "_15 output file" << endl;
        return 1;
    }
    out16.open(output_file_path + "_16", ofstream::trunc);
    if (not out16.is_open()) {
        cerr << "could not open " + output_file_path + "_16 output file" << endl;
        return 1;
    }
    out17.open(output_file_path + "_17", ofstream::trunc);
    if (not out17.is_open()) {
        cerr << "could not open " + output_file_path + "_17 output file" << endl;
        return 1;
    }
    out18.open(output_file_path + "_18", ofstream::trunc);
    if (not out18.is_open()) {
        cerr << "could not open " + output_file_path + "_18 output file" << endl;
        return 1;
    }
    out19.open(output_file_path + "_19", ofstream::trunc);
    if (not out19.is_open()) {
        cerr << "could not open " + output_file_path + "_19 output file" << endl;
        return 1;
    }
    out20.open(output_file_path + "_20", ofstream::trunc);
    if (not out20.is_open()) {
        cerr << "could not open " + output_file_path + "_20 output file" << endl;
        return 1;
    }
    out21.open(output_file_path + "_21", ofstream::trunc);
    if (not out21.is_open()) {
        cerr << "could not open " + output_file_path + "_21 output file" << endl;
        return 1;
    }
    out22.open(output_file_path + "_22", ofstream::trunc);
    if (not out22.is_open()) {
        cerr << "could not open " + output_file_path + "_22 output file" << endl;
        return 1;
    }
    out23.open(output_file_path + "_23", ofstream::trunc);
    if (not out23.is_open()) {
        cerr << "could not open " + output_file_path + "_23 output file" << endl;
        return 1;
    }
    out24.open(output_file_path + "_24", ofstream::trunc);
    if (not out24.is_open()) {
        cerr << "could not open " + output_file_path + "_24 output file" << endl;
        return 1;
    }
    out25.open(output_file_path + "_25", ofstream::trunc);
    if (not out25.is_open()) {
        cerr << "could not open " + output_file_path + "_25 output file" << endl;
        return 1;
    }
    out26.open(output_file_path + "_26", ofstream::trunc);
    if (not out26.is_open()) {
        cerr << "could not open " + output_file_path + "_26 output file" << endl;
        return 1;
    }
    out27.open(output_file_path + "_27", ofstream::trunc);
    if (not out27.is_open()) {
        cerr << "could not open " + output_file_path + "_27 output file" << endl;
        return 1;
    }
    out28.open(output_file_path + "_28", ofstream::trunc);
    if (not out28.is_open()) {
        cerr << "could not open " + output_file_path + "_28 output file" << endl;
        return 1;
    }
    out29.open(output_file_path + "_29", ofstream::trunc);
    if (not out29.is_open()) {
        cerr << "could not open " + output_file_path + "_29 output file" << endl;
        return 1;
    }
    out30.open(output_file_path + "_30", ofstream::trunc);
    if (not out30.is_open()) {
        cerr << "could not open " + output_file_path + "_30 output file" << endl;
        return 1;
    }
    out31.open(output_file_path + "_31", ofstream::trunc);
    if (not out31.is_open()) {
        cerr << "could not open " + output_file_path + "_31 output file" << endl;
        return 1;
    }
    out32.open(output_file_path + "_32", ofstream::trunc);
    if (not out32.is_open()) {
        cerr << "could not open " + output_file_path + "_32 output file" << endl;
        return 1;
    }
    out33.open(output_file_path + "_33", ofstream::trunc);
    if (not out33.is_open()) {
        cerr << "could not open " + output_file_path + "_33 output file" << endl;
        return 1;
    }
    out34.open(output_file_path + "_34", ofstream::trunc);
    if (not out34.is_open()) {
        cerr << "could not open " + output_file_path + "_34 output file" << endl;
        return 1;
    }
    out35.open(output_file_path + "_35", ofstream::trunc);
    if (not out35.is_open()) {
        cerr << "could not open " + output_file_path + "_35 output file" << endl;
        return 1;
    }
    out36.open(output_file_path + "_36", ofstream::trunc);
    if (not out36.is_open()) {
        cerr << "could not open " + output_file_path + "_36 output file" << endl;
        return 1;
    }
    out37.open(output_file_path + "_37", ofstream::trunc);
    if (not out37.is_open()) {
        cerr << "could not open " + output_file_path + "_37 output file" << endl;
        return 1;
    }
    out38.open(output_file_path + "_38", ofstream::trunc);
    if (not out38.is_open()) {
        cerr << "could not open " + output_file_path + "_38 output file" << endl;
        return 1;
    }
    out39.open(output_file_path + "_39", ofstream::trunc);
    if (not out39.is_open()) {
        cerr << "could not open " + output_file_path + "_39 output file" << endl;
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
        //marking the position in the editing channel as an edited adenosine
        else {
            dli_it->second.edit[pos_in_gene - 1] = '2';
        }
    }
    edit_in.close();
    //going over all the entries in the map
    int j = 0;
    for (map<string, dl_info>::iterator dli_it = deep_learning_info.begin(); dli_it != deep_learning_info.end(); dli_it++) {
        //going over all nucleotides in the sequence
        for (int i = 0; i < dli_it->second.edit.length(); i++) {
            //checking if the nuclaotide is a non-edited adenosine
            if (dli_it->second.edit[i] == '1') {
                //choosing which file the particular position is outputed to
                switch (j++ % 40) {
                    case 0:
                        out0 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 1:
                        out1 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 2:
                        out2 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 3:
                        out3 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 4:
                        out4 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 5:
                        out5 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 6:
                        out6 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 7:
                        out7 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 8:
                        out8 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 9:
                        out9 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 10:
                        out10 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 11:
                        out11 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 12:
                        out12 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 13:
                        out13 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 14:
                        out14 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 15:
                        out15 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 16:
                        out16 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 17:
                        out17 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 18:
                        out18 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 19:
                        out19 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 20:
                        out20 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 21:
                        out21 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 22:
                        out22 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 23:
                        out23 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 24:
                        out24 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 25:
                        out25 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 26:
                        out26 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 27:
                        out27 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 28:
                        out28 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 29:
                        out29 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 30:
                        out30 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 31:
                        out31 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 32:
                        out32 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 33:
                        out33 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 34:
                        out34 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 35:
                        out35 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 36:
                        out36 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 37:
                        out37 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 38:
                        out38 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                    case 39:
                        out39 << dli_it->first << "\t" << i + 1 << endl;
                        break;
                }
            }
        }
    }
    out0.close();
    out1.close();
    out2.close();
    out3.close();
    out4.close();
    out5.close();
    out6.close();
    out7.close();
    out8.close();
    out9.close();
    out10.close();
    out11.close();
    out12.close();
    out13.close();
    out14.close();
    out15.close();
    out16.close();
    out17.close();
    out18.close();
    out19.close();
    out20.close();
    out21.close();
    out22.close();
    out23.close();
    out24.close();
    out25.close();
    out26.close();
    out27.close();
    out28.close();
    out29.close();
    out30.close();
    out31.close();
    out32.close();
    out33.close();
    out34.close();
    out35.close();
    out36.close();
    out37.close();
    out38.close();
    out39.close();
}

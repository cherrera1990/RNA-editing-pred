/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a file with a list of positions within genes as given by the programs full_premrna_sequence.cpp and non_editing_adenosines.cpp and the secondary structure annotationn file as given by 
 * the structure_features.cpp program and gives a table where for each position there's a series of descriptors about the secondary structure in relation to the position. The table has a header giving all the
 * descriptors. Global descriptors give information about the structure of the whole gene, local descriptors give information within a window defined by an option and XClosest descriptors give detailed
 * about the N instances of each type closest to the annotated position (N defined by an option).
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

//a struct containing info for double strand fragments
typedef struct {
    int start1;
    int start2;
    int end1;
    int end2;
    int size;
} ds_info;

//a struct containing info for a particular feature of the secondary structure
typedef struct {
    char type;
    bool nick;
    int start1;
    int end1;
    int mid1;
    int size1;
    int start2;
    int end2;
    int mid2;
    int size2;
} feat_info;

//a struct containing info about a particular feature in relation to the target adenosine
typedef struct {
    char type;
    bool nick;
    bool inside;
    int closest_distance;
    int closest_rel_pos;
    int rel_pos_mid1;
    int rel_pos_mid2;
    int size1;
    int size2;
} local_feat_info;

//a struct containing info about a gene and all the global secondary structure
typedef struct {
    string gene_name;
    vector <ds_info> ds;
    vector <feat_info> feat;
    double dsPerc;
    double avgHLsize = 0;
    double avgILsize = 0;
    double avgBLsize = 0;
    int numHL = 0;
    int numIL = 0;
    int numBL = 0;
    int maxHLsize = 0;
    int maxILsize = 0;
    int maxBLsize = 0;
    int maxDS = 0;
} gene_info;

//a string comparer to use within maps
struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};

int main (int argc, char* argv[]) {
    //checking if the -h option is given
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << " RNA_editing_file RNA_structures_annotation_file [options]" << endl;
            cout << "Options: -s followed by the size of the sequence to be explored at both sides of each editing event (default = 50)   -m followed by the number of features of each type closest to the editing event to consider" << endl;
            return 0;
        }
    }
    //checking for the correct number of arguments
    if (argc < 3 or argc > 7) {
        cerr << "Usage: " << argv[0] << " RNA_editing_file RNA_structures_annotation_file" << endl;
        cerr << "Options: -s followed by the size of the sequence to be explored at both sides of each editing event (default = 50)   -m followed by the number of features of each type closest to the editing event to consider" << endl;
        return 1;
    }
    //implementing options
    int wsize = 50;
    int nfeat = 3;
    for (int i = 3; i < argc; i++) {
        string opt(argv[i]);
        if (opt.compare("-s") == 0 and i + 1 < argc) wsize = atoi(argv[i++ + 1]);
        else if (opt.compare("-m") == 0 and i + 1 < argc) nfeat = atoi(argv[i++ + 1]);
        else {
            cerr << "Usage: " << argv[0] << " RNA_editing_file RNA_structures_annotation_file" << endl;
            cerr << "Options: -s followed by the size of the sequence to be explored at both sides of each editing event (default = 50)   -m followed by the number of features of each type closest to the editing event to consider" << endl;
            return 1;
        }
    }
    //getting the input file paths from argumentsº
    const string edit_in_path(argv[1]);
    const string structures_in_path(argv[2]);
    //opening the structure features annotation file
    ifstream edit_in;
    ifstream structures_in;
    structures_in.open(structures_in_path.c_str());
    if (not structures_in.is_open()) {
        cerr << "could not open RNA structures annotation file" << endl;
        return 1;
    }
    //declaring the maps where we will store the info about each gene and the structure features for each gene
    string structures_line; 
    map<string, gene_info, comparer> genes_struct;
    map<string, gene_info>::iterator genstruct_it;
    //reading the structure features file line by line
    while (getline(structures_in, structures_line)) {
        //if the line is not a start of a new gene, then it is a feature annotation
        if (structures_line[0] != '>') {
            //we parse the line, check which type of feature is annotated and save the relevant info in the map
            vector<string> structures_parsed = split(structures_line, ' ');
            if (structures_parsed[0].compare("DS") == 0) {
                ds_info currDS;
                vector<string> p1 = split(structures_parsed[1], '=');
                vector<string> p1_parsed = split(p1[1], ',');
                currDS.start1 = atoi(p1_parsed[0].c_str());
                currDS.start2 = atoi(p1_parsed[1].c_str());
                vector<string> p2 = split(structures_parsed[2], '=');
                vector<string> p2_parsed = split(p2[1], ',');
                currDS.end1 = atoi(p2_parsed[0].c_str());
                currDS.end2 = atoi(p2_parsed[1].c_str());
                vector<string> s = split(structures_parsed[3], '=');
                currDS.size = atoi(s[1].c_str());
                genstruct_it->second.ds.push_back(currDS);
            }
            else if (structures_parsed[0].compare("SSMF") == 0 or structures_parsed[0].compare("SSBF") == 0 or structures_parsed[0].compare("SSEF") == 0) {
                feat_info currFeat;
                currFeat.type = 's';
                currFeat.nick = false;
                vector<string> st = split(structures_parsed[1], '=');
                currFeat.start1 = atoi(st[1].c_str());
                vector<string> en = split(structures_parsed[2], '=');
                currFeat.end1 = atoi(en[1].c_str());
                vector<string> s = split(structures_parsed[3], '=');
                currFeat.size1 = atoi(s[1].c_str());
                genstruct_it->second.feat.push_back(currFeat);
            }
            else if (structures_parsed[0].compare("HL") == 0) {
                feat_info currFeat;
                currFeat.type = 'h';
                currFeat.nick = false;
                vector<string> st = split(structures_parsed[1], '=');
                currFeat.start1 = atoi(st[1].c_str());
                vector<string> en = split(structures_parsed[2], '=');
                currFeat.end1 = atoi(en[1].c_str());
                vector<string> m = split(structures_parsed[3], '=');
                currFeat.mid1 = atoi(m[1].c_str());
                vector<string> s = split(structures_parsed[4], '=');
                currFeat.size1 = atoi(s[1].c_str());
                genstruct_it->second.feat.push_back(currFeat);
            }
            //in previous versions, we tried annotating branching points, but it had problems, now we annotate the single strand fragments between branches as SSMF
            /*else if (structures_parsed[0].compare("BP") == 0) {
                feat_info currFeat;
                currFeat.type = 'p';
                currFeat.nick = false;
                vector<string> st1 = split(structures_parsed[4], '=');
                currFeat.start1 = atoi(st1[1].c_str());
                vector<string> en1 = split(structures_parsed[5], '=');
                currFeat.end1 = atoi(en1[1].c_str());
                vector<string> s1 = split(structures_parsed[6], '=');
                currFeat.size1 = atoi(s1[1].c_str());
                vector<string> st2 = split(structures_parsed[7], '=');
                currFeat.start2 = atoi(st2[1].c_str());
                vector<string> en2 = split(structures_parsed[8], '=');
                currFeat.end2 = atoi(en2[1].c_str());
                vector<string> s2 = split(structures_parsed[9], '=');
                currFeat.size2 = atoi(s2[1].c_str());
                vector<string> st3 = split(structures_parsed[10], '=');
                currFeat.start3 = atoi(st3[1].c_str());
                vector<string> en3 = split(structures_parsed[11], '=');
                currFeat.end3 = atoi(en3[1].c_str());
                vector<string> s3 = split(structures_parsed[12], '=');
                currFeat.size3 = atoi(s3[1].c_str());
                genstruct_it->second.feat.push_back(currFeat);
            }*/
            else if (structures_parsed[0].compare("IL") == 0 or structures_parsed[0].compare("ILN") == 0) {
                feat_info currFeat;
                currFeat.type = 'i';
                currFeat.nick = (structures_parsed[0].compare("ILN") == 0);
                vector<string> st1 = split(structures_parsed[1], '=');
                currFeat.start1 = atoi(st1[1].c_str());
                vector<string> en1 = split(structures_parsed[2], '=');
                currFeat.end1 = atoi(en1[1].c_str());
                vector<string> m1 = split(structures_parsed[3], '=');
                currFeat.mid1 = atoi(m1[1].c_str());
                vector<string> s1 = split(structures_parsed[4], '=');
                currFeat.size1 = atoi(s1[1].c_str());
                vector<string> st2 = split(structures_parsed[5], '=');
                currFeat.start2 = atoi(st2[1].c_str());
                vector<string> en2 = split(structures_parsed[6], '=');
                currFeat.end2 = atoi(en2[1].c_str());
                vector<string> m2 = split(structures_parsed[7], '=');
                currFeat.mid2 = atoi(m2[1].c_str());
                vector<string> s2 = split(structures_parsed[8], '=');
                currFeat.size2 = atoi(s2[1].c_str());
                genstruct_it->second.feat.push_back(currFeat);
            }
            else if (structures_parsed[0].compare("BL") == 0 or structures_parsed[0].compare("BLN") == 0) {
                feat_info currFeat;
                currFeat.type = 'b';
                currFeat.nick = (structures_parsed[0].compare("BLN") == 0);
                vector<string> st1 = split(structures_parsed[1], '=');
                currFeat.start1 = atoi(st1[1].c_str());
                vector<string> en1 = split(structures_parsed[2], '=');
                currFeat.end1 = atoi(en1[1].c_str());
                vector<string> m1 = split(structures_parsed[3], '=');
                currFeat.mid1 = atoi(m1[1].c_str());
                vector<string> s1 = split(structures_parsed[4], '=');
                currFeat.size1 = atoi(s1[1].c_str());
                vector<string> st2 = split(structures_parsed[5], '=');
                currFeat.start2 = atoi(st2[1].c_str());
                vector<string> en2 = split(structures_parsed[6], '=');
                currFeat.end2 = atoi(en2[1].c_str());
                vector<string> m2 = split(structures_parsed[7], '=');
                currFeat.mid2 = atoi(m2[1].c_str());
                vector<string> s2 = split(structures_parsed[8], '=');
                currFeat.size2 = atoi(s2[1].c_str());
                genstruct_it->second.feat.push_back(currFeat);
            }
            else if (structures_parsed[0].compare("DSstats") == 0) {
                vector<string> p = split(structures_parsed[1], '=');
                genstruct_it->second.dsPerc = atof(p[1].c_str());
                vector<string> max = split(structures_parsed[2], '=');
                genstruct_it->second.maxDS = atoi(max[1].c_str());
            }
            else if (structures_parsed[0].compare("HLstats") == 0) {
                vector<string> num = split(structures_parsed[1], '=');
                genstruct_it->second.numHL = atoi(num[1].c_str());
                vector<string> av = split(structures_parsed[2], '=');
                genstruct_it->second.avgHLsize = atof(av[1].c_str());
                vector<string> max = split(structures_parsed[3], '=');
                genstruct_it->second.maxHLsize = atoi(max[1].c_str());
            }
            else if (structures_parsed[0].compare("ILstats") == 0) {
                vector<string> num = split(structures_parsed[1], '=');
                genstruct_it->second.numIL = atoi(num[1].c_str());
                vector<string> av = split(structures_parsed[2], '=');
                genstruct_it->second.avgILsize = atof(av[1].c_str());
                vector<string> max = split(structures_parsed[3], '=');
                genstruct_it->second.maxILsize = atoi(max[1].c_str());
            }
            else if (structures_parsed[0].compare("BLstats") == 0) {
                vector<string> num = split(structures_parsed[1], '=');
                genstruct_it->second.numBL = atoi(num[1].c_str());
                vector<string> av = split(structures_parsed[2], '=');
                genstruct_it->second.avgBLsize = atof(av[1].c_str());
                vector<string> max = split(structures_parsed[3], '=');
                genstruct_it->second.maxBLsize = atoi(max[1].c_str());
            }
            else cerr << "Warning: line \"" << structures_line << "\" not recognised" << endl; 
        }
        else {
            //if the line corresponds to a new gene, we create a new map entry for the gene
            structures_line.erase(0, 1);
            vector<string> structures_parsed = split(structures_line, '\t');
            gene_info curr_gene_info;
            curr_gene_info.gene_name = structures_parsed[0];
            curr_gene_info.ds = vector<ds_info> ();
            curr_gene_info.feat = vector<feat_info> ();
            genstruct_it = genes_struct.insert(pair<string, gene_info> (structures_parsed[1], curr_gene_info)).first;
        }
    }
    structures_in.close();
    //opening the editing input file
    edit_in.open(edit_in_path.c_str());
    if (not edit_in.is_open()) {
        cerr << "could not open RNA editing file" << endl;
        return 1;
    }
    string line;
    //writing the first line as the header of the table
    cout << "TranscriptId\tPosInTranscript\tGlobalPercDS\tGlobalMaxDSSize\tGlobalNumHL\tGlobalMaxHLSize\tGlobalAvgHLSize\tGlobalNumIL\tGlobalMaxILSize\tGlobalAvgILSize\tGlobalNumBL\tGlobalMaxBLSize\tGlobalAvgBLSize\tEventInDS\tEventInNick\tLocalDSperc\tLocalDistToDS\tLocalClosestDSSize\tLocalNumHL\tLocalMaxHLSize\tLocalAvgHLSize\tLocalNumIL\tLocalMaxILSize\tLocalAvgILSize\tLocalNumBL\tLocalMaxBLSize\tLocalAvgBLSize\tlocalNumSS\tLocalMaxSSSize\tLocalAvgSSSize";
    for (int i = 1; i <= nfeat; i++) {
        cout << "\t" << i << "ClosestHLEventInside\t" << i << "ClosestHLDistanceToEvent\t" << i << "ClosestHLClosestEndRelPos\t" << i << "ClosestHLMidRelPos\t" << i << "ClosestHLSize";
    }
    for (int i = 1; i <= nfeat; i++) {
        cout << "\t" << i << "ClosestILEventInside\t" << i << "ClosestILIsNick\t" << i << "ClosestILDistanceToEvent\t" << i << "ClosestILClosestEventStrandEndRelPos\t" << i << "ClosestILEventStrandMidRelPos\t" << i << "ClosestILEventStrandSize\t" << i << "ClosestILOppositeStrandMidRelPos\t" << i << "ClosestILOppositeStrandSize";
    }
    for (int i = 1; i <= nfeat; i++) {
        cout << "\t" << i << "ClosestBLEventInside\t" << i << "ClosestBLIsNick\t" << i << "ClosestBLDistanceToEvent\t" << i << "ClosestBLClosestEventStrandEndRelPos\t" << i << "ClosestBLEventStrandMidRelPos\t" << i << "ClosestBLEventStrandSize\t" << i << "ClosestBLOppositeStrandMidRelPos\t" << i << "ClosestBLOppositeStrandSize";
    }
    for (int i = 1; i <= nfeat; i++) {
        cout << "\t" << i << "ClosestSSEventInside\t" << i << "ClosestSSDistanceToEvent\t" << i << "ClosestSSClosestEndRelPos\t" << i << "ClosestSSSize";
    }
    cout << endl;
    //reading the editing input file line by line
    while (getline(edit_in, line)) {
        //parsing the line and getting the position of the target adenosine
        vector<string> parsed = split(line, '\t');
        string gene_id = parsed[1];
        int pos_in_gene = atoi(parsed[2].c_str());
        //finding the entry for the relevant gene
        map<string, gene_info>::iterator genstruct_it = genes_struct.find(gene_id);
        if (genstruct_it == genes_struct.end()) {
            cerr << "Warning: Something went wrong, Unable to find gene " << gene_id << " in structures file" << endl; 
        }
        else {
            //calculating the start and end of the window around the target adenosine 
            int wstart = pos_in_gene - wsize;
            int wend = pos_in_gene + wsize;
            //declaring and initialising variables relevant to the DS fragments
            bool in_DS = false;
            int dist_to_DS = -1;
            int closest_DS_length = -1;
            int sum_DS_length = 0;
            bool in_nick = false;
            //first, we go through all the double strand fragments in the gene
            for (int i = 0; i < genstruct_it->second.ds.size(); i++) {
                //we check, for both "strands", how much of the DS fragment is inside the window and add the corresponding length to the total DS length inside the window
                if (genstruct_it->second.ds[i].start1 >= wstart and genstruct_it->second.ds[i].end1 <= wend) sum_DS_length += genstruct_it->second.ds[i].end1 - genstruct_it->second.ds[i].start1 + 1;
                else if (genstruct_it->second.ds[i].start1 <= wstart and genstruct_it->second.ds[i].end1 >= wend) sum_DS_length += wend - wstart + 1;
                else if (genstruct_it->second.ds[i].start1 <= wstart and genstruct_it->second.ds[i].end1 >= wstart) sum_DS_length += genstruct_it->second.ds[i].end1 - wstart + 1;
                else if (genstruct_it->second.ds[i].start1 <= wend and genstruct_it->second.ds[i].end1 >= wend) sum_DS_length += wend - genstruct_it->second.ds[i].start1 + 1;
                if (genstruct_it->second.ds[i].end2 >= wstart and genstruct_it->second.ds[i].start2 <= wend) sum_DS_length += genstruct_it->second.ds[i].start2 - genstruct_it->second.ds[i].end2 + 1;
                else if (genstruct_it->second.ds[i].end2 <= wstart and genstruct_it->second.ds[i].start2 >= wend) sum_DS_length += wend - wstart + 1;
                else if (genstruct_it->second.ds[i].end2 <= wstart and genstruct_it->second.ds[i].start2 >= wstart) sum_DS_length += genstruct_it->second.ds[i].start2 - wstart + 1;
                else if (genstruct_it->second.ds[i].end2 <= wend and genstruct_it->second.ds[i].start2 >= wend) sum_DS_length += wend - genstruct_it->second.ds[i].end2 + 1;
                //we check if the target adenosine is inside the current DS fragment and add the corresponding information
                if (pos_in_gene >= genstruct_it->second.ds[i].start1 and pos_in_gene <= genstruct_it->second.ds[i].end1 or pos_in_gene >= genstruct_it->second.ds[i].end2 and pos_in_gene <= genstruct_it->second.ds[i].start2) {
                    in_DS = true;
                    dist_to_DS = 0;
                    closest_DS_length = genstruct_it->second.ds[i].end1 - genstruct_it->second.ds[i].start1 + 1;
                }
                //if the target adenosine is not inside the current fragment but is inside the window,  we calculate the smallest distance from the DS fragment to the target adenosine and if it is the closest distance yet, we update the closest distance to DS
                else if (genstruct_it->second.ds[i].end1 >= wstart and genstruct_it->second.ds[i].end1 <= wend or genstruct_it->second.ds[i].start1 >= wstart and genstruct_it->second.ds[i].start1 <= wend
                    or genstruct_it->second.ds[i].end2 >= wstart and genstruct_it->second.ds[i].end2 <= wend or genstruct_it->second.ds[i].start2 >= wstart and genstruct_it->second.ds[i].start2 <= wend) {
                    int smallest_dist_to_pos = abs(genstruct_it->second.ds[i].start1 - pos_in_gene);
                    int dist_to_pos = abs(genstruct_it->second.ds[i].end1 - pos_in_gene);
                    if (dist_to_pos < smallest_dist_to_pos) smallest_dist_to_pos = dist_to_pos;
                    dist_to_pos = abs(genstruct_it->second.ds[i].start2 - pos_in_gene);
                    if (dist_to_pos < smallest_dist_to_pos) smallest_dist_to_pos = dist_to_pos;
                    dist_to_pos = abs(genstruct_it->second.ds[i].end2 - pos_in_gene);
                    if (dist_to_pos < smallest_dist_to_pos) smallest_dist_to_pos = dist_to_pos;
                    if (smallest_dist_to_pos < dist_to_DS or dist_to_DS == -1) {
                        dist_to_DS = smallest_dist_to_pos;
                        closest_DS_length = genstruct_it->second.ds[i].end1 - genstruct_it->second.ds[i].start1 + 1;
                    }
                    
                    /*if (genstruct_it->second.ds[i].end1 >= wstart and genstruct_it->second.ds[i].end1 <= wend and genstruct_it->second.ds[i].start1 < wstart) sum_DS_length += genstruct_it->second.ds[i].end1 - wstart + 1;
                    else if (genstruct_it->second.ds[i].start1 >= wstart and genstruct_it->second.ds[i].start1 <= wend and genstruct_it->second.ds[i].end1 > wend) sum_DS_length += wend - genstruct_it->second.ds[i].start1 + 1;    
                    else if (genstruct_it->second.ds[i].start2 >= wstart and genstruct_it->second.ds[i].start2 <= wend and genstruct_it->second.ds[i].end2 < wstart) sum_DS_length += genstruct_it->second.ds[i].start2 - wstart + 1;
                    else if (genstruct_it->second.ds[i].end2 >= wstart and genstruct_it->second.ds[i].end2 <= wend and genstruct_it->second.ds[i].start2 > wend) sum_DS_length += wend - genstruct_it->second.ds[i].end2 + 1;
                    if (genstruct_it->second.ds[i].start1 >= wstart and genstruct_it->second.ds[i].end1 <= wend)  sum_DS_length += genstruct_it->second.ds[i].end1 - genstruct_it->second.ds[i].start1 + 1;
                    if (genstruct_it->second.ds[i].end2 >= wstart and genstruct_it->second.ds[i].start2 <= wend)  sum_DS_length += genstruct_it->second.ds[i].start2 - genstruct_it->second.ds[i].end2 + 1;*/
                }
            }
            //declaring and initialising variables relevant for all the features that are not DS
            double localDSperc = (double) sum_DS_length / (double) (wend - wstart + 1);
            int localNumHL = 0;
            int localNumIL = 0;
            int localNumBL = 0;
            int localNumSS = 0;
            vector<local_feat_info> localfeats = vector<local_feat_info> ();
            vector<local_feat_info> closest_HL = vector<local_feat_info> ();
            vector<local_feat_info> closest_IL = vector<local_feat_info> ();
            vector<local_feat_info> closest_BL = vector<local_feat_info> ();
            vector<local_feat_info> closest_SS = vector<local_feat_info> ();
            //we go over each feature, chekc its type and calculate the relevant information, saving it in the vectors
            for (int i = 0; i < genstruct_it->second.feat.size(); i++) {
                if (genstruct_it->second.feat[i].type == 's') {
                    local_feat_info localf;
                    if (pos_in_gene >= genstruct_it->second.feat[i].start1 and pos_in_gene <= genstruct_it->second.feat[i].end1) {
                        localf.type = 's';
                        localf.inside = true;
                        localf.nick = false;
                        localf.closest_distance = 0;
                        int dist1 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist2 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        if (abs(dist1) < abs(dist2)) localf.closest_rel_pos = dist1;
                        else localf.closest_rel_pos = dist2;
                        localf.size1 = genstruct_it->second.feat[i].size1;
                        localfeats.push_back(localf);
                        localNumSS++;
                    }
                    else {
                        localf.type = 's';
                        localf.inside = false;
                        localf.nick = false;
                        int dist1 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist2 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        if (abs(dist1) < abs(dist2)) localf.closest_rel_pos = dist1;
                        else localf.closest_rel_pos = dist2;
                        localf.closest_distance = abs(localf.closest_rel_pos);
                        localf.size1 = genstruct_it->second.feat[i].size1;
                        if (genstruct_it->second.feat[i].start1 >= wstart and genstruct_it->second.feat[i].start1 <= wend or genstruct_it->second.feat[i].end1 >= wstart and genstruct_it->second.feat[i].end1 <= wend) {
                            localfeats.push_back(localf);
                            localNumSS++;
                        }
                    }
                    int j = 0;
                    bool found = false;
                    while (j < closest_SS.size() and not found) {
                        if (closest_SS[j].closest_distance > localf.closest_distance) {
                            closest_SS.insert(closest_SS.begin() + j, localf);
                            found = true;
                        }
                        else j++;
                    }
                    if (closest_SS.size() < nfeat and not found) closest_SS.push_back(localf);
                    else if (closest_SS.size() > nfeat) closest_SS.pop_back();
                }
                else if (genstruct_it->second.feat[i].type == 'h') {
                    local_feat_info localf;
                    if (pos_in_gene >= genstruct_it->second.feat[i].start1 and pos_in_gene <= genstruct_it->second.feat[i].end1) {
                        localf.type = 'h';
                        localf.inside = true;
                        localf.nick = false;
                        localf.closest_distance = 0;
                        int dist1 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist2 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        if (abs(dist1) < abs(dist2)) localf.closest_rel_pos = dist1;
                        else localf.closest_rel_pos = dist2;
                        localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                        localf.size1 = genstruct_it->second.feat[i].size1;
                        localfeats.push_back(localf);
                        localNumHL++;
                    }
                    else {
                        localf.type = 'h';
                        localf.inside = false;
                        localf.nick = false;
                        int dist1 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist2 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        if (abs(dist1) < abs(dist2)) localf.closest_rel_pos = dist1;
                        else localf.closest_rel_pos = dist2;
                        localf.closest_distance = abs(localf.closest_rel_pos);
                        localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                        localf.size1 = genstruct_it->second.feat[i].size1;
                        if (genstruct_it->second.feat[i].start1 >= wstart and genstruct_it->second.feat[i].start1 <= wend or genstruct_it->second.feat[i].end1 >= wstart and genstruct_it->second.feat[i].end1 <= wend) {
                            localfeats.push_back(localf);
                            localNumHL++;
                        }
                    }
                    int j = 0;
                    bool found = false;
                    while (j < closest_HL.size() and not found) {
                        if (closest_HL[j].closest_distance > localf.closest_distance) {
                            closest_HL.insert(closest_HL.begin() + j, localf);
                            found = true;
                        }
                        else j++;
                    }
                    if (closest_HL.size() < nfeat and not found) closest_HL.push_back(localf);
                    else if (closest_HL.size() > nfeat) closest_HL.pop_back();
                }
                else if (genstruct_it->second.feat[i].type == 'i') {
                    local_feat_info localf;
                    if (pos_in_gene >= genstruct_it->second.feat[i].start1 and pos_in_gene <= genstruct_it->second.feat[i].end1 or pos_in_gene >= genstruct_it->second.feat[i].end2 and pos_in_gene <= genstruct_it->second.feat[i].start2) {
                        localf.type = 'i';
                        localf.inside = true;
                        localf.nick = genstruct_it->second.feat[i].nick;
                        if (localf.nick) in_nick = true;
                        localf.closest_distance = 0;
                        int dist11 =  genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist12 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        int dist1;
                        if (abs(dist11) < abs(dist12)) dist1 = dist11;
                        else dist1 = dist12;
                        int dist21 = genstruct_it->second.feat[i].start2 - pos_in_gene;
                        int dist22 = genstruct_it->second.feat[i].end2 - pos_in_gene;
                        int dist2;
                        if (abs(dist21) < abs(dist22)) dist2 = dist21;
                        else dist2 = dist22;
                        if (abs(dist1) < abs(dist2)) {
                            localf.closest_rel_pos = dist1;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size1;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size2;
                        }
                        else {
                            localf.closest_rel_pos = dist2;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size2;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size1;
                        }
                        localfeats.push_back(localf);
                        localNumIL++;
                    }
                    else  {
                        localf.type = 'i';
                        localf.inside = false;
                        localf.nick = genstruct_it->second.feat[i].nick;
                        int dist11 =  genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist12 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        int dist1;
                        if (abs(dist11) < abs(dist12)) dist1 = dist11;
                        else dist1 = dist12;
                        int dist21 = genstruct_it->second.feat[i].start2 - pos_in_gene;
                        int dist22 = genstruct_it->second.feat[i].end2 - pos_in_gene;
                        int dist2;
                        if (abs(dist21) < abs(dist22)) dist2 = dist21;
                        else dist2 = dist22;
                        if (abs(dist1) < abs(dist2)) {
                            localf.closest_rel_pos = dist1;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size1;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size2;
                        }
                        else {
                            localf.closest_rel_pos = dist2;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size2;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size1;
                        }
                        localf.closest_distance = abs(localf.closest_rel_pos);
                        if (genstruct_it->second.feat[i].start1 >= wstart and genstruct_it->second.feat[i].start1 <= wend or genstruct_it->second.feat[i].end1 >= wstart and genstruct_it->second.feat[i].end1 <= wend
                            or genstruct_it->second.feat[i].start2 >= wstart and genstruct_it->second.feat[i].start2 <= wend or genstruct_it->second.feat[i].end2 >= wstart and genstruct_it->second.feat[i].end2 <= wend) {
                            localfeats.push_back(localf);
                            localNumIL++;
                        }
                    }
                    int j = 0;
                    bool found = false;
                    while (j < closest_IL.size() and not found) {
                        if (closest_IL[j].closest_distance > localf.closest_distance) {
                            closest_IL.insert(closest_IL.begin() + j, localf);
                            found = true;
                        }
                        else j++;
                    }
                    if (closest_IL.size() < nfeat and not found) closest_IL.push_back(localf);
                    else if (closest_IL.size() > nfeat) closest_IL.pop_back();
                }
                else if (genstruct_it->second.feat[i].type == 'b') {
                    local_feat_info localf;
                    if (pos_in_gene >= genstruct_it->second.feat[i].start1 and pos_in_gene <= genstruct_it->second.feat[i].end1 or pos_in_gene >= genstruct_it->second.feat[i].end2 and pos_in_gene <= genstruct_it->second.feat[i].start2) {
                        localf.type = 'b';
                        localf.inside = true;
                        localf.nick = genstruct_it->second.feat[i].nick;
                        if (localf.nick) in_nick = true;
                        localf.closest_distance = 0;
                        int dist11 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist12 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        int dist1;
                        if (abs(dist11) < abs(dist12)) dist1 = dist11;
                        else dist1 = dist12;
                        int dist21 = genstruct_it->second.feat[i].start2 - pos_in_gene;
                        int dist22 = genstruct_it->second.feat[i].end2 - pos_in_gene;
                        int dist2;
                        if (abs(dist21) < abs(dist22)) dist2 = dist21;
                        else dist2 = dist22;
                        if (abs(dist1) < abs(dist2)) {
                            localf.closest_rel_pos = dist1;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size1;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size2;
                        }
                        else {
                            localf.closest_rel_pos = dist2;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size2;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size1;
                        }
                        localfeats.push_back(localf);
                        localNumBL++;
                    }
                    else {
                        localf.type = 'b';
                        localf.inside = false;
                        localf.nick = genstruct_it->second.feat[i].nick;
                        int dist11 = genstruct_it->second.feat[i].start1 - pos_in_gene;
                        int dist12 = genstruct_it->second.feat[i].end1 - pos_in_gene;
                        int dist1;
                        if (abs(dist11) < abs(dist12)) dist1 = dist11;
                        else dist1 = dist12;
                        int dist21 = genstruct_it->second.feat[i].start2 - pos_in_gene;
                        int dist22 = genstruct_it->second.feat[i].end2 - pos_in_gene;
                        int dist2;
                        if (abs(dist21) < abs(dist22)) dist2 = dist21;
                        else dist2 = dist22;
                        if (abs(dist1) < abs(dist2)) {
                            localf.closest_rel_pos = dist1;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size1;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size2;
                        }
                        else {
                            localf.closest_rel_pos = dist2;
                            localf.rel_pos_mid1 = genstruct_it->second.feat[i].mid2 - pos_in_gene;
                            localf.size1 = genstruct_it->second.feat[i].size2;
                            localf.rel_pos_mid2 = genstruct_it->second.feat[i].mid1 - pos_in_gene;
                            localf.size2 = genstruct_it->second.feat[i].size1;
                        }
                        localf.closest_distance = abs(localf.closest_rel_pos);
                        if (genstruct_it->second.feat[i].start1 >= wstart and genstruct_it->second.feat[i].start1 <= wend or genstruct_it->second.feat[i].end1 >= wstart and genstruct_it->second.feat[i].end1 <= wend
                            or genstruct_it->second.feat[i].start2 >= wstart and genstruct_it->second.feat[i].start2 <= wend or genstruct_it->second.feat[i].end2 >= wstart and genstruct_it->second.feat[i].end2 <= wend) {
                            localfeats.push_back(localf);
                            localNumBL++;
                        }
                    }
                    int j = 0;
                    bool found = false;
                    while (j < closest_BL.size() and not found) {
                        if (closest_BL[j].closest_distance > localf.closest_distance) {
                            closest_BL.insert(closest_BL.begin() + j, localf);
                            found = true;
                        }
                        else j++;
                    }
                    if (closest_BL.size() < nfeat and not found) closest_BL.push_back(localf);
                    else if (closest_BL.size() > nfeat) closest_BL.pop_back();
                }
            }
            //we now calculate the relevant stats for the local window
            int sumLocalHLSize = 0;
            int localHLMaxSize = -1;
            int sumLocalILSize = 0;
            int localILMaxSize = -1;
            int sumLocalBLSize = 0;
            int localBLMaxSize = -1;
            int sumLocalSSSize = 0;
            int localSSMaxSize = -1;
            for (int i = 0; i < localfeats.size(); i++) {
                if (localfeats[i].type == 's') {
                    if (localfeats[i].size1 > localSSMaxSize) localSSMaxSize = localfeats[i].size1;
                    sumLocalSSSize += localfeats[i].size1;
                }
                if (localfeats[i].type == 'h') {
                    if (localfeats[i].size1 > localHLMaxSize) localHLMaxSize = localfeats[i].size1;
                    sumLocalHLSize += localfeats[i].size1;
                }
                if (localfeats[i].type == 'i') {
                    if (localfeats[i].size1 > localILMaxSize) localILMaxSize = localfeats[i].size1;
                    sumLocalILSize += localfeats[i].size1;
                    if (localfeats[i].size2 > localILMaxSize) localILMaxSize = localfeats[i].size2;
                    sumLocalILSize += localfeats[i].size2;
                }
                if (localfeats[i].type == 'b') {
                    if (localfeats[i].size1 > localBLMaxSize) localBLMaxSize = localfeats[i].size1;
                    sumLocalILSize += localfeats[i].size1;
                    if (localfeats[i].size2 > localBLMaxSize) localBLMaxSize = localfeats[i].size2;
                    sumLocalBLSize += localfeats[i].size2;
                }
            }
            double localHLAvgSize, localILAvgSize, localBLAvgSize, localSSAvgSize;
            /*if (localNumHL == 0) localHLAvgSize = -1;
            else*/ localHLAvgSize = (double) sumLocalHLSize / (double) localNumHL;
            /*if (localNumIL == 0) localILAvgSize = -1;
            else*/ localILAvgSize = (double) sumLocalILSize / (double) (localNumIL * 2);
            /*if (localNumBL == 0) localBLAvgSize = -1;
            else*/ localBLAvgSize = (double) sumLocalBLSize / (double) (localNumBL * 2);
            /*if (localNumSS == 0) localSSAvgSize = -1;
            else*/ localSSAvgSize = (double) sumLocalSSSize / (double) localNumSS;
            //we output all the descriptors for the current target adenosine
            cout << gene_id << "\t" << pos_in_gene << "\t" << genstruct_it->second.dsPerc << "\t" << genstruct_it->second.maxDS << "\t" << genstruct_it->second.numHL << "\t" << genstruct_it->second.maxHLsize  << "\t" << genstruct_it->second.avgHLsize << "\t" << genstruct_it->second.numIL << "\t" << genstruct_it->second.maxILsize  << "\t" << genstruct_it->second.avgILsize << "\t" << genstruct_it->second.numBL << "\t" << genstruct_it->second.maxBLsize  << "\t" << genstruct_it->second.avgBLsize << "\t" << in_DS << "\t" << in_nick << "\t" << localDSperc << "\t" << dist_to_DS  << "\t" << closest_DS_length << "\t" << localNumHL << "\t" << localHLMaxSize << "\t" << localHLAvgSize << "\t" << localNumIL << "\t" << localILMaxSize << "\t" << localILAvgSize << "\t" << localNumBL << "\t" << localBLMaxSize << "\t" << localBLAvgSize << "\t" << localNumSS << "\t" << localSSMaxSize << "\t" << localSSAvgSize;
            for (int i = 1; i <= nfeat; i++) {
                if (i - 1 < closest_HL.size()) cout << "\t" << closest_HL[i - 1].inside << "\t" << closest_HL[i - 1].closest_distance << "\t" << closest_HL[i - 1].closest_rel_pos << "\t" << closest_HL[i - 1].rel_pos_mid1 << "\t" << closest_HL[i - 1].size1;
                else cout << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1;
            }
            for (int i = 1; i <= nfeat; i++) {
                if (i - 1 < closest_IL.size()) cout << "\t" << closest_IL[i - 1].inside << "\t" << closest_IL[i - 1].nick << "\t" << closest_IL[i - 1].closest_distance << "\t" << closest_IL[i - 1].closest_rel_pos << "\t" << closest_IL[i - 1].rel_pos_mid1 << "\t" << closest_IL[i - 1].size1 << "\t" << closest_IL[i - 1].rel_pos_mid2 << "\t" << closest_IL[i - 1].size2;
                else cout << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1;
            }
            for (int i = 1; i <= nfeat; i++) {
                if (i - 1 < closest_BL.size()) cout << "\t" << closest_BL[i - 1].inside << "\t" << closest_BL[i - 1].nick << "\t" << closest_BL[i - 1].closest_distance << "\t" << closest_BL[i - 1].closest_rel_pos << "\t" << closest_BL[i - 1].rel_pos_mid1 << "\t" << closest_BL[i - 1].size1 << "\t" << closest_BL[i - 1].rel_pos_mid2 << "\t" << closest_BL[i - 1].size2;
                else cout << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1;
            }
            for (int i = 1; i <= nfeat; i++) {
                if (i - 1 < closest_SS.size()) cout << "\t" << closest_SS[i - 1].inside << "\t" << closest_SS[i - 1].closest_distance << "\t" << closest_SS[i - 1].closest_rel_pos << "\t" << closest_SS[i - 1].size1;
                else cout << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1;
            }
            cout << endl;
        }
    }
    edit_in.close();
}

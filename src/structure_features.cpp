/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program takes as input a file with secondary structures in a dot-parenthesis format as given by`programs such as linearfold and gives as an output an annotation of the features of each secondary structure.
 * The program detects the following features: Hairpin loop (HL), Inner loop (IL), Bulge loop (BL), Double strand fragments (DS) and non-loop single strand fragments (SSBF, SSMF, SSEF). The output has the following
 * format: for each sequence, first a line with the title as given in the input file (must begin with '>' using the fasta format convention), followed by information about each feature, one feature per line. Generally
 * for each feature, the start and end pair is given and the size of the feature (or two sizes in the case of BL and IL). For BL and IL, there are specific annotations for features of size 1, which we call nicks and
 * are annotated as ILN and BLN. After all the features of a gene are listed, there are a few lines giving a few statistics for that gene.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stack>

using namespace std;

int main (int argc, char* argv[]) {
    //checking for the -h option
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "structures_file" << endl;
            return 0;
        }
    }
    //checking the correct number of arguments
    if (argc != 2) {
            cerr << "Usage: " << argv[0] << "structures_file" << endl;
            return 1;
    }
    //getting the input file path from the argument
    const string str_in_path(argv[1]);
    //opening the input file
    ifstream str_in;
    str_in.open(str_in_path.c_str());
    if (not str_in.is_open()) {
        cerr << "could not open RNA structures file" << endl;
        return 1;
    }
    //we use this vector to store the genes we have encountered, in order to rename them if they have already appeared
    vector<string> gene_ids;
    //we declare the variables that will hold the three lines for each gene
    string title;
    string seq;
    string str;
    //for each gene, we read the title, the sequence and the structure prediction
    while (getline(str_in, title)) {
        if (not getline(str_in, seq)) {
            cerr << "error " << title << " has no sequence" << endl;
            return 1;
        }
        if (not getline(str_in, str)) {
            cerr << "error " << title << " has no structure" << endl;
            return 1;
        }
        //the structure line has at the end a score separated by a space, we remove it
        int ws_pos = str.find_first_of(' ');
        str.erase(ws_pos);
        //we create a vector that will hold the positions of all the pairs of base-paired nucleotides
        vector<pair<int, int> > dsPosAnnot;
        //for efficiency purposes, we reserve memory for this vector
        dsPosAnnot.reserve(str.length());
        //saving the paired parenthesis positions
        int lparPos = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str[i] == '(') {
                dsPosAnnot.push_back(pair<int, int> (i, -1));
                lparPos++;
            }
            else if (str[i] == ')') {
                int j = lparPos - 1;
                while (dsPosAnnot[j].second != -1) j--;
                dsPosAnnot[j].second = i;
            }
        }
        //if we find the current gene already annotated, we change the current name adding _altN to it
        //note that in the full_premrna_sequence program we already add the chromosome or scaffold to the gene name
        //so this step should be redundant
        string orig_title = title;
        int alts = 0;
        for (int i = 0; i < gene_ids.size(); i++) {
            if (gene_ids[i].compare(title) == 0) title = orig_title + "_alt" + to_string(++alts);
        }
        cout << title << endl;
        gene_ids.push_back(title);
        //we calculate the percentage of nucleotides in DS
        double dsPerc = (double) (lparPos * 2) / (double) str.length();
        //we declare all the variables that will hold general statistics for this sequence's structure
        int sumHLsize = 0;
        int sumILsize = 0;
        int sumBLsize = 0;
        int numHL = 0;
        int numIL = 0;
        int numBL = 0;
        int maxHLsize = 0;
        int maxILsize = 0;
        int maxBLsize = 0;
        int maxDS = 0;
        //if we detect no DSRNA, then there are no features in the structure
        if (dsPerc == 0.0) cout << "NO DSRNA" << endl;
        else {
            //we use this stack to keep track of the double strand branches or hairpins we encounter: a branch starts with the end of a SSBF or SSMF and ends with a hairpin loop (HL).
            stack<pair<int, int> > branch_starts;
            //we get the start of the and end first DS fragment
            branch_starts.push(dsPosAnnot[0]);
            int currDSStart = dsPosAnnot[0].first;
            int currDSEnd = dsPosAnnot[0].second;
            //anything before the start of the first DS fragment is a SSBF
            cout << "SSBF start=1 end=" << currDSStart << " size=" << currDSStart << endl; //"\t";
            //we keep track of the start of the current double strand fragment and the corresponding base-paired position
            int currDSLStart = dsPosAnnot[0].first;
            int currDSRstart = dsPosAnnot[0].second;
            //we go over all the base-paired positions
            for (int i = 1; i < dsPosAnnot.size(); i++) {
                //we compare the current pair, to the next pair
                int currLPar = dsPosAnnot[i - 1].first;
                int nextLPar = dsPosAnnot[i].first;
                int currRPar = dsPosAnnot[i - 1].second;
                int nextRPar = dsPosAnnot[i].second;
                //if the current pair is not strictly contiguous to the next pair on both strands, we have a break in the DS and a feature
                if (currLPar + 1 != nextLPar or currRPar - 1 != nextRPar) {
                    //we calculate and annotate the DS fragment that has just ended
                    int DSsize = currLPar - currDSLStart + 1;
                    if (DSsize > maxDS) maxDS = DSsize;
                    cout << "DS startPair=" << currDSLStart + 1 << "," << currDSRstart + 1 << " endPair=" << currLPar + 1 << "," << currRPar + 1 << " size=" << DSsize << endl; //"\t";
                    //we update the start of the next DS fragmet
                    currDSLStart = nextLPar;
                    currDSRstart = nextRPar;
                    //in order to recognise the type of structure, we need to check if the next DS fragment is a continuation of the current one, by looking for the next DS pair that could have continuity
                    int greatestRPar = -1;
                    int greatestLPar = -1;
                    for (int j = i; j < dsPosAnnot.size(); j++) {
                        if (dsPosAnnot[j].second < currRPar and dsPosAnnot[j].second > greatestRPar) {
                            greatestRPar = dsPosAnnot[j].second;
                            greatestLPar = dsPosAnnot[j].first;
                        }
                    }
                    //if there's no continuity in the DS, but the next open parenthesis continues has the continuation of the current right parenthesis, it is a hairpin loop
                    if (greatestRPar == -1 and nextLPar > currRPar) {
                        //we calculate and print all the relevant stats for the hairpin loop
                        int start = currLPar + 1;
                        int end = currRPar - 1;
                        int size = end - start + 1;
                        int mid = start + (size/2);
                        numHL++;
                        sumHLsize += size;
                        if (size > maxHLsize) maxHLsize = size;
                        cout << "HL srart=" << start + 1 << " end=" << end + 1 << " mid=" << mid + 1 << " size=" << size << endl; //"\t";
                        //we check how many branches are closed by this hairpin loop and annotate the corresponding SSMF
                        pair<int, int> last_branch_start = branch_starts.top();
                        branch_starts.pop();
                        pair<int, int> current_branch_start;
                        if (not branch_starts.empty()) current_branch_start = branch_starts.top();
                        pair<int, int> next_branch_start = dsPosAnnot[i];
                        while (next_branch_start.first > current_branch_start.second and not branch_starts.empty()) {
                            int j = 0;
                            while (last_branch_start.second != dsPosAnnot[j].second) j++;
                            int min = str.length();
                            for (int k = 0; k < j; k++) {
                                if (dsPosAnnot[k].second < min and dsPosAnnot[k].second > last_branch_start.second) {
                                    min = dsPosAnnot[k].second;
                                }
                            }
                            int ssmf_start = last_branch_start.second + 1;
                            int ssmf_end = min - 1;
                            int ssmf_size = ssmf_end - ssmf_start + 1;
                            cout << "SSMF start=" << ssmf_start + 1 << " end=" << ssmf_end + 1 << " size=" << ssmf_size << endl;
                            branch_starts.pop();
                            last_branch_start = current_branch_start;
                            if (not branch_starts.empty()) current_branch_start = branch_starts.top();
                        }
                        int ssmf_start = last_branch_start.second + 1;
                        int ssmf_end = next_branch_start.first - 1;
                        int ssmf_size = ssmf_end - ssmf_start + 1;
                        cout << "SSMF start=" << ssmf_start + 1 << " end=" << ssmf_end + 1 << " size=" << ssmf_size << endl;
                        //we push the start of the next branch
                        branch_starts.push(next_branch_start);
                    }
                    //if there's strict continuity in the DS fragments, then it is either an inner loop or a bulge loop
                    else if (greatestRPar == nextRPar) {
                        //we calculate all the stats for the current feature, bulge loops and inner loops have the same stats
                        int start1 = currLPar + 1;
                        int end1 = nextLPar - 1;
                        int size1 = end1 - start1 + 1;
                        int mid1 = start1 + (size1/2);
                        int start2 = currRPar - 1;
                        int end2 = nextRPar + 1;
                        int size2 = start2 - end2 + 1;
                        int mid2 = end2 + (size2/2);
                        //if one of the two sizes is 0, then the feature is a bulge loop (BL), and if the other size is 1, then it is a bulge loop nick (BLN)
                        if (size1 == 0 or size2 == 0) {
                            if (size1 == 1 or size2 == 1) {
                                cout << "BLN srart1=" << start1 + 1 << " end1=" << end1 + 1 << " mid1=" << mid1 + 1 << " size1=" << size1 << " srart2=" << start2 + 1 << " end2=" << end2 + 1 << " mid2=" << mid2 + 1 << " size2=" << size2 << endl; //"\t";
                            }
                            else {
                                cout << "BL srart1=" << start1 + 1 << " end1=" << end1 + 1 << " mid1=" << mid1 + 1 << " size1=" << size1 << " srart2=" << start2 + 1 << " end2=" << end2 + 1 << " mid2=" << mid2 + 1 << " size2=" << size2 << endl; //"\t";
                            }
                            numBL++;
                            sumBLsize += size1 + size2;
                            if (size1 + size2 > maxBLsize) maxBLsize = size1 + size2;
                        }
                        //if both sizes are 1, it is an inner loop nick (ILN)
                        else if (size1 == 1 and size2 == 1) {
                            cout << "ILN srart1=" << start1 + 1 << " end1=" << end1 + 1 << " mid1=" << mid1 + 1 << " size1=" << size1 << " srart2=" << start2 + 1 << " end2=" << end2 + 1 << " mid2=" << mid2 + 1 << " size2=" << size2 << endl; //"\t";
                            numIL++;
                            sumILsize += size1 + size2;
                            if (size1 + size2 > maxILsize) maxILsize = size1 + size2;
                        }
                        //if all sizes are greater than 0 and one is greater than one, it is an inner loop
                        else {
                            cout << "IL srart1=" << start1 + 1 << " end1=" << end1 + 1 << " mid1=" << mid1 + 1 << " size1=" << size1 << " srart2=" << start2 + 1 << " end2=" << end2 + 1 << " mid2=" << mid2 + 1 << " size2=" << size2 << endl; //"\t";
                            numIL++;
                            sumILsize += size1 + size2;
                            if (size1 + size2 > maxILsize) maxILsize = size1 + size2;
                        }
                    }
                    //if the DS has no direct continuity and it is not a HL then it is a SSMF
                    else if (greatestRPar != -1) {
                        branch_starts.push(dsPosAnnot[i]);
                        int start1 = currLPar + 1;
                        int end1 = nextLPar - 1;
                        int size1 = end1 - start1 + 1;
                        //cout << "BP pairing1=" << currLPar + 1 << "," << currRPar + 1 << " pairing2=" << nextLPar +  1 << "," << nextRPar + 1 << " pairing3=" << greatestLPar + 1 << "," << greatestRPar + 1 << endl; //"\t";
                        cout << "SSMF start=" << start1 + 1 << " end=" << end1 + 1 << " size=" << size1 << endl;
                    }
                    //in case we find a structure that doesn't fit any of the previous features, we annotate it as UNKNOWN_FEATURE
                    else {
                        cout << "UNKNOWN_FEATURE pairing1=" << currLPar + 1 << "," << currRPar + 1 << " pairing2=" << nextLPar + 1 << "," << nextRPar + 1 << endl; //"\t";
                    }
                }
            }
            //after analysing all the structure, we calculate and output all the statistics for the current sequence
            int DSsize = dsPosAnnot[dsPosAnnot.size() - 1].first - currDSLStart + 1;
            if (DSsize > maxDS) maxDS = DSsize;
            cout << "DS startPair=" << currDSLStart + 1 << "," << currDSRstart + 1 << " endPair=" << dsPosAnnot[dsPosAnnot.size() - 1].first + 1 << "," << dsPosAnnot[dsPosAnnot.size() - 1].second + 1 << " size=" << DSsize << endl; //"\t";
            int start = dsPosAnnot[dsPosAnnot.size() - 1].first + 1;
            int end = dsPosAnnot[dsPosAnnot.size() - 1].second - 1;
            int size = end - start + 1;
            int mid = start + (size/2);
            numHL++;
            sumHLsize += size;
            if (size > maxHLsize) maxHLsize = size;
            cout << "HL srart=" << start + 1 << " end=" << end + 1 << " mid=" << mid + 1 << " size=" << size << endl; //"\t";
            pair<int, int> last_branch_start = branch_starts.top();
            branch_starts.pop();
            pair<int, int> current_branch_start;
            if (not branch_starts.empty())  current_branch_start = branch_starts.top();
            while (not branch_starts.empty()) {
                int j = 0;
                while (last_branch_start.second != dsPosAnnot[j].second) j++;
                int min = str.length();
                for (int k = 0; k < j; k++) {
                    if (dsPosAnnot[k].second < min and dsPosAnnot[k].second > last_branch_start.second) {
                        min = dsPosAnnot[k].second;
                    }
                }
                int ssmf_start = last_branch_start.second + 1;
                int ssmf_end = min - 1;
                int ssmf_size = ssmf_end - ssmf_start + 1;
                cout << "SSMF start=" << ssmf_start + 1 << " end=" << ssmf_end + 1 << " size=" << ssmf_size << endl;
                branch_starts.pop();
                last_branch_start = current_branch_start;
                if (not branch_starts.empty()) current_branch_start = branch_starts.top();
            }
            cout << "SSEF start=" << last_branch_start.second + 2 << " end=" << str.length() << " size=" << str.length() - 1 - last_branch_start.second << endl; //"\t";
            double avgHLsize = 0;
            double avgILsize = 0;
            double avgBLsize = 0;
            if (numHL != 0) avgHLsize = (double) sumHLsize / (double) numHL;
            if (numIL != 0) avgILsize = (double) sumILsize / (double) numIL;
            if (numBL != 0) avgBLsize = (double) sumBLsize / (double) numBL;
            cout << "DSstats perc=" << dsPerc << " maxSize=" << maxDS << endl; //"\t"
            cout << "HLstats num=" << numHL << " avgsize=" << avgHLsize << " maxHLsize=" << maxHLsize << endl; //"\t"
            cout << "ILstats num=" << numIL << " avgsize=" << avgILsize << " maxILsize=" << maxILsize << endl; //"\t"
            cout << "BLstats num=" << numBL << " avgsize=" << avgBLsize << " maxBLsize=" << maxBLsize << endl; //"\t"
            //cout << endl;
        }
    }
}
/*
 * Author Michal Zawisza Alvarez, contact m.zawisza@ub.edu
 * This program uses as input a list of (edited) genomic positions, a reference genome and a gene annotation flle (either a gtf or gff) to output all the premrna sequences
 * that contain the edited positions as well as the same list of positions, but adding the gene identifier nad position within the premrna ssequence. The sequences and new
 * posistion lines are printed in separate files, and there's an option to divide the sequecnes among a number of files. If the gene annotation file have a gene entry within
 * the whole locus annotated, there's an option to use the gene notation directly. Otherwise, it uses the transcript entries to get the minimum starting and maximum ending
 * positions of each gene to extract the premrna sequence. If the gene is annotated in the negative strand, this program outputs the complementary sequence, but does not
 * do the reverse. There are options to ignore sequences over a given length or with a percentage of gaps over a given treshold.
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

// gtf_info is a struct used to store all the relevant information of an annotated gene, including the premrna sequence
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

//a functions that transforms a sequence into its reverse complementary
string reverse_complementary_seq(const string &seq) {
    string cmp(seq);
    char c;
    for (int i = 0; i < seq.length(); i++) {
        c = seq[i];
        if (c == 'A') cmp[cmp.size() - i - 1] = 'T';
        else if (c == 'a') cmp[cmp.size() - i - 1] = 't';
        else if (c == 'T') cmp[cmp.size() - i - 1] = 'A';
        else if (c == 't') cmp[cmp.size() - i - 1] = 'a';
        else if (c == 'G') cmp[cmp.size() - i - 1] = 'C';
        else if (c == 'g') cmp[cmp.size() - i - 1] = 'c';
        else if (c == 'C') cmp[cmp.size() - i - 1] = 'G';
        else if (c == 'c') cmp[cmp.size() - i - 1] = 'g';
    }
    return cmp;
}

//the main function of the program
int main (int argc, char* argv[]) {
    //checking if the -h option is given
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << " filtered_variants_size genome_fasta_file gtf_file sequence_output_file editing_output_file [options]" << endl;
            cout << "Options: -d followed by the number of files to divide the results    -s followed by the transcript size treshold over which the genes will be ignored    -n followed by the percentage of Ns in a transcript over which a gene will be ignored     -c followed by the maximum chromosome or scaffold size (for efficiency reasons)     -g if entire gene locus are annotated in the gtf or gff (one gene entry per gene id)    -f if the gene annotation file is in gff format (gtf by default    -p if the fasta header contains extra information after the chromosome/scaffold id" << endl;
            return 0;
        }
    }
    //checking minimum number of arguments
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " filtered_variants_size genome_fasta_file gtf_file sequence_output_file editing_output_file [options]" << endl;
        cerr << "Options: -d followed by the number of files to divide the results    -s followed by the transcript size treshold over which the genes will be ignored    -n followed by the percentage of Ns in a transcript over which a gene will be ignored     -c followed by the maximum chromosome or scaffold size (for efficiency reasons)     -g if entire gene locus are annotated in the gtf or gff (one gene entry per gene id)    -f if the gene annotation file is in gff format (gtf by default    -p if the fasta header contains extra information after the chromosome/scaffold id" << endl;
        return 1;
    }
    //getting the input and output file paths from the first 5 arguments
    const string var_in_path(argv[1]);
    const string fasta_in_path(argv[2]);
    const string gtf_in_path(argv[3]);
    const string seq_out_path(argv[4]);
    const string edit_out_path(argv[5]);
    //initialising optional variables with default values
    int div = -1;
    int size_treshold = -1;
    int n_perc = -1;
    int max_chr = 536870912;
    bool annotated_genes = false;
    bool gff = false;
    bool parse_fasta_title = false;
    //checking for valid options and getting the corresponding values
    for (int i = 6; i < argc; i++) {
        string opt(argv[i]);
        if (i < argc) {
            if (opt.compare("-d") == 0) div = atoi(argv[++i]);
            else if (opt.compare("-s") == 0) size_treshold = atoi(argv[++i]);
            else if (opt.compare("-n") == 0) n_perc = atoi(argv[++i]);
            else if (opt.compare("-c") == 0) max_chr = atoi(argv[++i]);
            else if (opt.compare("-g") == 0) annotated_genes = true;
            else if (opt.compare("-f") == 0) gff = true;
            else if (opt.compare("-p") == 0) parse_fasta_title = true;
            else {
                cerr << "Usage: " << argv[0] << " filtered_variants_size genome_fasta_file gtf_file sequence_output_file editing_output_file [options]" << endl;
                cerr << "Options: -d followed by the number of files to divide the results    -s followed by the transcript size treshold over which the genes will be ignored    -n followed by the percentage of Ns in a transcript over which a gene will be ignored     -c followed by the maximum chromosome or scaffold size (for efficiency reasons)     -g if entire gene locus are annotated in the gtf or gff (one gene entry per gene id)    -f if the gene annotation file is in gff format (gtf by default    -p if the fasta header contains extra information after the chromosome/scaffold id" << endl;
                return 1;
            }
        }
        else {
            cerr << "Usage: " << argv[0] << " filtered_variants_size genome_fasta_file gtf_file sequence_output_file editing_output_file [options]" << endl;
            cerr << "Options: -d followed by the number of files to divide the results    -s followed by the transcript size treshold over which the genes will be ignored    -n followed by the percentage of Ns in a transcript over which a gene will be ignored     -c followed by the maximum chromosome or scaffold size (for efficiency reasons)     -g if entire gene locus are annotated in the gtf or gff (one gene entry per gene id)    -f if the gene annotation file is in gff format (gtf by default    -p if the fasta header contains extra information after the chromosome/scaffold id" << endl;
            return 1;
        }
    }
    //opening the positions list file and the genome fasta file
    ifstream var_in;
    ifstream fasta_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open filtered variants file" << endl;
        return 1;
    }
    fasta_in.open(fasta_in_path.c_str());
    if (not fasta_in.is_open()) {
        cerr << "could not open fasta file" << endl;
        return 1;
    }
    //variables to read and store the genome file, we store the genome sequences in a map with each chromosome or scaffold as the key
    string fasta_line; 
    map<string, string, comparer> genome_seq;
    map<string, string>::iterator genseq_it;
    bool eof = false;
    bool end_scaffold = false;
    //reading the genome fasta file line by line
    while (getline(fasta_in, fasta_line)) {
        //if the line is not a title, we read the line and append it to the current sequence
        if (fasta_line[0] != '>') {
            genseq_it->second += fasta_line;
        }
        //if the line is a title, we create the new entry in the map and reserve the space given as the length of the largest chromosome for efficiency purposes
        else {
            
            fasta_line.erase(0, 1);
            if (parse_fasta_title) {
                vector<string> fasta_parsed = split(fasta_line, ' ');
                fasta_line = fasta_parsed[0];
            }
            genseq_it = genome_seq.insert(pair<string, string> (fasta_line, string())).first;
            genseq_it->second.reserve(max_chr);
        }
    }
    fasta_in.close();
    //we open the gtf or gff file
    ifstream gtf_in;
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open transcript reference file" << endl;
        return 1;
    }
    //we will save the relevant gene annotation info in a map of a vector of the corresponding struct, using the scaffold or chromosome as key
    map<string, vector<gtf_info>, comparer> scaf_gtf_info;
    string gtf_line;
    //we read the gtf or gff file line by line
    while (getline(gtf_in, gtf_line)) {
        // we ignore all the lines that begin with a #
        if (gtf_line[0] != '#') {
            //we parse the gtf line using tabs
            vector<string> gtf_parsed = split(gtf_line, '\t');
            // we check if the line corresponds to a gene notation if the -g option is active or a transcript notation otherwise
            bool is_gene;
            if (annotated_genes) is_gene = (gtf_parsed[2].compare("gene") == 0);
            else is_gene = (gtf_parsed[2].compare("transcript") == 0);
            if (is_gene) {
                //we get the chromosome or scaffold id from the parsed line
                string scaf_id = gtf_parsed[0];
                //we create the struct which will contain the relevant information about this entry
                gtf_info n_info;
                //we parse the last field into subfields using ;
                vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ';');
                vector<string> gtf_subfield_parsed;
                //we look for the subfield that has the gene name
                for (int i = 0; i < gtf_last_field_parsed.size(); i++) {
                    //gtf and gff use sligthly different formats for the subfields
                    if (gff) gtf_subfield_parsed = split(gtf_last_field_parsed[i], '=');
                    else {
                        if (gtf_last_field_parsed[i][0] == ' ') gtf_last_field_parsed[i].erase(0,1);
                        gtf_subfield_parsed = split(gtf_last_field_parsed[i], ' ');
                    }
                    if (gtf_subfield_parsed[0].compare("gene_id") == 0 or gtf_subfield_parsed[0].compare("ID") == 0) n_info.gene_id = gtf_subfield_parsed[1];
                }
                if (not gff) {
                    n_info.gene_id.erase(n_info.gene_id.length() - 1, 1);
                    n_info.gene_id.erase(0, 1);
                }
                //we use the gene_id as gene_name for information purposes
                n_info.gene_name = n_info.gene_id;
                //but as identifier we use the gene_id and the chromosome or scaffold in conjuction to avoid problems with alternative notations
                n_info.gene_id += "_" + scaf_id;
                //we get the annotated start, end and the strand information from the annotation line
                n_info.gtf_start = atoi(gtf_parsed[3].c_str());
                n_info.gtf_end = atoi(gtf_parsed[4].c_str());
                n_info.positive_strand = (gtf_parsed[6][0] == '+');
                //the insert function of map inserts a new entry if ihat enttry is not in the map, but if one already exists with the given key, it doesn't insert but returns the existing entry
                map<string, vector<gtf_info> >::iterator gtf_it = scaf_gtf_info.insert(pair<string, vector<gtf_info> > (scaf_id, vector<gtf_info> ())).first;
                //if the -g option is not used we check if the entry of the annotated gene exists
                if (not annotated_genes) {
                    bool found = false;
                    int i = 0;
                    while (not found and i < gtf_it->second.size()) {
                        if (n_info.gene_name.compare(gtf_it->second[i].gene_name) == 0) found = true;
                        else i++;
                    }
                    //if the entry for the gene exists, we update the start and end position if necessary
                    if (found) {
                        if (n_info.gtf_start < gtf_it->second[i].gtf_start) gtf_it->second[i].gtf_start = n_info.gtf_start;
                        if (n_info.gtf_end < gtf_it->second[i].gtf_end) gtf_it->second[i].gtf_end = n_info.gtf_end;
                    }
                    //if he entry does not exist, we add the entry as is
                    else {
                        gtf_it->second.push_back(n_info);
                    }
                }
                //if we use the -g option, we assume there's only one notation for each chromosome or scaffold and add the new entry directly
                else gtf_it->second.push_back(n_info);
            }
        }
    }
    gtf_in.close();
    string line;
    //this map will contain the info of all the genes that contain at least one editing target, including the premnra sequence and all the lines that correspond to the editing targets within the gene
    map<string, gtf_info, comparer> gtf_and_seqs;
    //reading the variants or RNA editing annotation file line by line
    while (getline(var_in, line)) {
        //parsing the line using tabs
        vector<string> parsed = split(line, '\t');
        //checking if the right changes are annotated for A-to-I editing
        if (parsed[2][0] == 'T' and parsed[3][0] == 'C' and parsed[4][0] == '-' or parsed[2][0] == 'A' and parsed[3][0] == 'G' and parsed[4][0] == '+') {
            //retrieving the information about thegenomic sequence and the genes annotated in the scaffold 
            string scaf_id = parsed[0];
            map<string, string>::iterator fas_it = genome_seq.find(scaf_id);
            if (fas_it == genome_seq.end()) {
                cerr << "Something went wrong, Unable to find scaffold " << scaf_id << " in fasta file" << endl; 
                return 1;
            }
            map<string, vector<gtf_info> >::iterator gtf_it = scaf_gtf_info.find(parsed[0]);
            if (gtf_it == scaf_gtf_info.end()) {
                cerr << "Warning: Something went wrong, Unable to find scaffold " << parsed[0] << " in gtf file" << endl; 
            }
            else {
                //parsing the position and strand of the current editing annotation
                int pos = atoi(parsed[1].c_str());
                bool positive_strand = (parsed[4][0] == '+');
                //looking for the gene that contains the current position
                bool found = false;
                for (int i = 0; i < gtf_it->second.size(); i++) {
                    int start = gtf_it->second[i].gtf_start;
                    int end = gtf_it->second[i].gtf_end;
                    //when we find it, we take the whole sequence from the genome
                    if (pos >= start and pos <= end and positive_strand == gtf_it->second[i].positive_strand) {
                        if (not found) found = true;
                        string seq(fas_it->second, start - 1, end - start + 1);
                        //if the target is annotated in the negative strand, then we calculate the complementary sequence (but not the reverse)
                        if (not positive_strand) seq = reverse_complementary_seq(seq);
                        //we calculate the local position of the target inside the sequence
                        int trans_pos = pos - start + 1;
                        //we search for the current gene in the map
                        map<string, gtf_info>::iterator gtf_seqs_it = gtf_and_seqs.find(gtf_it->second[i].gene_id);
                        //if the gene hasn't been saved yet in the map, we create the corresponding entry and save it
                        if (gtf_seqs_it == gtf_and_seqs.end()) {
                            gtf_info n_info;
                            n_info.gene_id = gtf_it->second[i].gene_id;
                            n_info.gene_name = gtf_it->second[i].gene_name;
                            n_info.gtf_start = start;
                            n_info.gtf_end = end;
                            n_info.positive_strand = gtf_it->second[i].positive_strand;
                            n_info.editing_lines = vector<string> ();
                            string seq(fas_it->second, start - 1, end - start + 1);
                            if (not positive_strand) seq = reverse_complementary_seq(seq);
                            n_info.seq = seq;
                            gtf_seqs_it = gtf_and_seqs.insert(pair<string, gtf_info> (n_info.gene_id, n_info)).first;
                        }
                        //we add the current editing target annotation to the list of targets within that gene
                        gtf_seqs_it->second.editing_lines.push_back(line);
                    }
                }
                if (not found) cerr << "Warning: Something went wrong, Unable to find a transcript that contains \"" << line << "\" in gtf file" << endl; 
            }
        }
    }
    var_in.close();
    //we go over the saved genes and if the appropriate options are active, we discard those genes that are too large or have too many Ns in the sequence
    map<string, gtf_info>::iterator gtf_seqs_it = gtf_and_seqs.begin();
    ofstream seq_out;
    ofstream edit_out;
    if (n_perc != -1 or size_treshold != -1) {
        int n_ignored = 0;
        int s_ignored = 0;
        while (gtf_seqs_it != gtf_and_seqs.end()) {
            map<string, gtf_info>::iterator next_it = next(gtf_seqs_it);
            bool ignore_n = false;
            bool ignore_s = false;
            if (n_perc != -1) {
                int n_count = 0;
                for (int i = 0; i < gtf_seqs_it->second.seq.length(); i++) {
                    if (gtf_seqs_it->second.seq[i] == 'N' or gtf_seqs_it->second.seq[i] == 'n') n_count++;
                }
                if ((double) n_count / (double) gtf_seqs_it->second.seq.length() * 100 > (double) n_perc) {
                    ignore_n = true;
                    n_ignored++;
                }
            }
            if (size_treshold != -1) {
                if (gtf_seqs_it->second.seq.length() > size_treshold) {
                    ignore_s = true;
                    s_ignored++;
                }
            }
            if (ignore_n or ignore_s) {
                gtf_and_seqs.erase(gtf_seqs_it);
            }
            gtf_seqs_it = next_it;
        }
        if (n_perc != -1) cerr << "genes ignored by N percentage: " << n_ignored << endl;
        if (size_treshold != -1) cerr << "genes ignored by size: " << s_ignored << endl;
    }
    gtf_seqs_it = gtf_and_seqs.begin();
    //if the corresponding option is active, we calculate the number of sequences in each of the separate files
    int divsize;
    if (div == -1) {
        seq_out.open(seq_out_path, ofstream::trunc);
        if (not seq_out.is_open()) {
            cerr << "could not open sequence output file" << endl;
            return 1;
        }
    }
    else divsize = gtf_and_seqs.size()/div;
    edit_out.open(edit_out_path, ofstream::trunc);
    if (not edit_out.is_open()) {
        cerr << "could not open editing output file" << endl;
        return 1;
    }
    //we now go over every gene left in the map, and print the sequence to one (or multiple) file and all the target lines with the local position added at the beginning, corresponding to targets within that gene to another file
    int i = 0;
    while (gtf_seqs_it != gtf_and_seqs.end()) {
        if (div != -1) {
            if (i%divsize == 0) {
                if (i != 0) seq_out.close();
                seq_out.open(seq_out_path + "_" + to_string(i/divsize + 1) + ".fa", ofstream::trunc);
                if (not seq_out.is_open()) {
                    cerr << "could not open sequence output file" << endl;
                    return 1;
                }
            }
        }
        for (int i = 0; i < gtf_seqs_it->second.editing_lines.size(); i++) {
            line = gtf_seqs_it->second.editing_lines[i];
            vector<string> parsed = split(line, '\t');
            int pos = atoi(parsed[1].c_str());
            if (gtf_seqs_it->second.positive_strand) trans_pos = pos - gtf_seqs_it->second.gtf_start + 1;
            else trans_pos = gtf_seqs_it->second.gtf_end - pos + 1;
            edit_out << gtf_seqs_it->second.gene_name << '\t' << gtf_seqs_it->second.gene_id << '\t' << trans_pos << '\t' << line << endl;
        }
        seq_out << ">" << gtf_seqs_it->second.gene_name << '\t' << gtf_seqs_it->second.gene_id << endl;
        seq_out << gtf_seqs_it->second.seq << endl;
        gtf_seqs_it++;
        i++;
    }
    seq_out.close();
    edit_out.close();
}


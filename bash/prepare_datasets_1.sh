#PART 1 of 6: please, make sure every process ends before running the next part of this script. When using the small dataset example, take into account that steps 4 to 6 will take
#a disproportionate amount of time only to add more examples to the negative dataset pool.
#The first argument of this script is the table with RNA-editing annotated positions. We provide a small dataset, but the full table for human and mouse can be downloaded from
#REDIPORTAL: http://srv00.recas.ba.infn.it/atlas/download.html. The second argument should be the genome fasta file (we use the hg38 human genome assembly for the small dataset example)
#Note that if the fasta file has more annotations after the chromosome or scaffold name, the -p option should be added to the full_premrna_sequence and full_premrna_not_edited_sequence 
#commands. The third argumnt should be the gtf or gff gene annotation file (we used the ncbiRefSeq gtf file for the human hg38 genome assembly for the small dataset example). Note that
#if a gff file is used, the -g and -f #options should be added to the full_premrna_sequence and full_premrna_not_edited_sequence commands and the locus_lengths_gff program should be used 
#instead of locus_lenghts. The fourth argument should be the prefix for all the files created by this script. The fifth argument should be the path for the linearfold executable. 
#Note that the small dataset we provide is a subset of the human REDIPORTAL dataset, and it will not achieve good results when using it for training and testing the Random Forest and Neural
#Networks models.
#For generating the table from the RNA-seq and DNA-seq data, you can use vcfs and our 2 auxiliary programs (variants_filter_gff_v2.cpp and DNA_reads_editing_filter.cpp). However, make sure
#you transform the vcf lines to match the format of the 6 first columns of the REDIPORTAL tables.
#After finishing this script, for Neural Networks you should use the ClipsAndCodeGenes.py python script (provided in the folder biLSTM/tools) with the _DL output files of this script and then
#notebooks for GoogleColab (provided in the folder biLSTM/notebooks). For Random Forest, you should use the R scripts (provided in the folder RF) with the _editing_structures_DB output files
#of this script.

#compiling all the necessary C++ programs. This only needs to be done once. If you want to run this script again, we recommend that you comment the compilation lines for a slightly faster performance.
g++ ../src/locus_lengths.cpp -o locus_lengths.x
g++ ../src/locus_lengths_gff.cpp -o locus_lengths_gff.x
g++ ../src/full_premrna_sequence.cpp -o full_premrna_sequence.x
g++ ../src/structure_features.cpp -o structure_features.x
g++ ../src/deep_learning_encoder.cpp -o deep_learning_encoder.x
g++ ../src/non_editing_adenosines.cpp -o non_editing_adenosines.x
g++ ../src/editing_structures_DB_v2.cpp -o editing_structures_DB_v2.x
g++ ../src/random_lines.cpp -o random_lines.x


#Getting the list of all genes annotated in the gtf
./locus_lengths.x $3 > $4_lengths
#After loading the lengths in R and using the summary option on the table, we can calculate the statistical outliers with the formula: Q3 + 1,5 * IQR and use the result
#when running the next program.
mkdir seqs
#we get the pre-mRNA sequences that contain the editing events and we add the local position of each event in the sequence to the annotation. We ignore all the sequences
#that have moore than 20% Ns, and that are longer than 50616 nucleotides.
./full_premrna_sequence.x $1 $2 $3 seqs/$4_premrna_seqs $4_editing_lines -d 60 -n 20 -s 50616 2> premrna_seqs.log
#the following script runs linearfold in parallel with the 60 sequence files to get the predicted secondary structure.
mkdir structures
bash bash_linearfold.sh $4 $5


#it is necessary to make sure for that all the runs of the previous script end before continuing running the next steps. Please, make sure all the linarfold processes end before running the next part.

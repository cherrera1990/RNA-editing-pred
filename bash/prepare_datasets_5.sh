#PART 5 of 6: please, make sure every process ends before running the next part of this script. When using the small dataset example, take into account that steps 4 to 6 will take
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

#it is necessary to make sure for that all the runs of the previous script end before continuing running the next steps. Please, make sure all the linarfold processes end before running the next part.
#the following script runs linearfold in parallel with the 60 sequence files to get the predicted secondary structure.
bash bash_regroup_files_NES.sh $4
#this program annotates the features of the predicted structures
./structure_features.x structures/$4_premrna_NES_struct_all > structures/$4_premrna_NES_struct_annot 2> features_NES.log
#generatiing an empy file to indicate that there are no edited adenosines in the sequences
touch empty_file
#this program generetaes the initial multi-channel encoding for the Neural Network models.
./deep_learning_encoder.x empty_file structures/$4_premrna_NES_struct_all structures/$4_premrna_NES_truct_annot > $4_NES_DL 2> $4_NES_DL.log
#this program generates the local positions of the adenosines that are not marked as edited in the initial dataset. This will be the adenosines marked as negative for the Random Forest datasets,
./non_editing_adenosines.x empty_file structures/$4_premrna_NES_struct_all $4_non_editing_NES_lines 2> $4_NES_non_editing.log
#these scripts generate the table output with the features of non_edited datasets that are used in the R scripts for the Random Forest models.
bash bash_not_edited_seqs_structures_DB_m5.sh $4

#it is necessary to make sure for that all the runs of the previous script end before continuing running the next steps. Please, make sure all the editing_structures_DB_v2 processes end before running
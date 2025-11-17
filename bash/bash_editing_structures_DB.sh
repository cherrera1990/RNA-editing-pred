nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 50 > $1_editing_structures_DB_v2_50_3 2> $1_DB_50_3.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 200 > $1_editing_structures_DB_v2_200_3 2> $1_DB_200_3.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 500 > $1_editing_structures_DB_v2_500_3 2> $1_DB_500_3.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 1000 > $1_editing_structures_DB_v2_1000_3 2> $1_DB_1000_3.log &

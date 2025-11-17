nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 50 -m 1 > $1_editing_structures_DB_v2_50_1 2> $1_DB_50_1.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 200 -m 1 > $1_editing_structures_DB_v2_200_1 2> $1_DB_200_1.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 500 -m 1 > $1_editing_structures_DB_v2_500_1 2> $1_DB_500_1.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 1000 -m 1 > $1_editing_structures_DB_v2_1000_1 2> $1_DB_1000_1.log &

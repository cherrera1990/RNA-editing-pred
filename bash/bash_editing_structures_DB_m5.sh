nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 50 -m 5 > $1_editing_structures_DB_v2_50_5 2> $1_DB_50_5.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 200 -m 5 > $1_editing_structures_DB_v2_200_5 2> $1_DB_200_5.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 500 -m 5 > $1_editing_structures_DB_v2_500_5 2> $1_DB_500_5.log &
nohup ./editing_structures_DB_v2.x $1_editing_lines structures/$1_premrna_struct_annot -s 1000 -m 5 > $1_editing_structures_DB_v2_1000_5 2> $1_DB_1000_5.log &

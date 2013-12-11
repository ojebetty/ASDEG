#!/bin/bash

##### Setting parameters #####
gene_file=$1
data_path=$2
trimming_path=$3
read_length=$4
case=$5
control=$6
script_path=$7
output=$8
mutithread=$9
########snp unit file########################

snp_chromsome=($(less $gene_file | cut -f 1))
snp_gene=($(less $gene_file | cut -f 3))
snp_coordinate_max=($(less $gene_file | cut -f 2))
snp_candidate_count=`less $gene_file | cut -f 2 | wc -l`

for (( i=0; i<$snp_candidate_count; i=i+1))

  do
     less $gene_file | grep -w "${snp_coordinate_max[i]}" > $trimming_path/${snp_gene[0]}_folder/snp_${snp_gene[0]}_${snp_coordinate_max[i]}

done

touch $trimming_path/log



for x in $trimming_path/${snp_gene[0]}_folder/snp_${snp_gene[0]}_*

 do

     snp_coordinate_max=`less $x | cut -f 2`
     gene_file=$x
     snp_chromsome=`less $x | cut -f 1`
     snp_gene=`less $x | cut -f 3`
    while [ $(ps -Af | grep "Trimming_function" | wc -l) -gt $mutithread ]

      do
         sleep 5
  done


      echo "data_path:$data_path"
      echo "trimming_path:$trimming_path"
      echo "script_path:$script_path"
      echo "read_length:$read_length"
      echo "case:$case"
      echo "snp_chromsome:$snp_chromsome"
      echo "snp_gene:$snp_gene"
      echo "snp_coordinate_max:$snp_coordinate_max"
      echo "gene_file:$gene_file"
      echo "output:$output"   
      ${script_path}/Trimming_function.sh $data_path $trimming_path $script_path $read_length $case $snp_chromsome $snp_gene $snp_coordinate_max $gene_file $output&
      ${script_path}/Trimming_function.sh $data_path $trimming_path $script_path $read_length $control $snp_chromsome $snp_gene $snp_coordinate_max $gene_file $output&

     sleep 2
done

sleep 1

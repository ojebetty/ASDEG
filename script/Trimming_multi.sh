#!/bin/bash

trimming_path=$2
data_path=$1
mutithread=$3
read_length=$4
case_file=$5
control_file=$6
script_path=$7
output=$8

gene_name=($(ls $trimming_path | grep "gene_name" | cut -d '_' -f 3 ))
total_gene=`ls $trimming_path | grep "gene_name" | cut -d '_' -f 3 | wc -l`

for ((i=0; i<$total_gene; i=i+1))

 do
    if ! [ -d $trimming_path/"${gene_name[i]}"_folder ]; then
 
       mkdir $trimming_path/"${gene_name[i]}"_folder
   fi

done

##### 利用基因的名稱作單位，開多顆cpu同時進行工作 #####

for x in $trimming_path/gene_name_*

 do
    while [ $(ps -Af | grep "Trimming_multi_snp" | wc -l) -gt $mutithread ]
       
      do
         sleep 5
  done

     ${script_path}/Trimming_multi_snp.sh $x $data_path $trimming_path $read_length $case_file $control_file $script_path $output $mutithread&
     sleep 1

done

sleep 5

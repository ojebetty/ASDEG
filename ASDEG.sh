#!/bin/bash

data_path=
script_path=
snp2gene=
thread=1
output_dir=ASDEG

while getopts p:S:G:t:c:C:g:o:h OPTION
do
    case $OPTION in
     p) data_path="$OPTARG" ;;
     S) script_path="$OPTARG"
        main=1 ;;
     G) snp2gene="$OPTARG" ;;    
     g) gene_name_column="$OPTARG" ;;
     t) thread="$OPTARG" ;;
     c) case="$OPTARG" ;;
     C) control="$OPTARG" ;;
     o) output_dir="$OPTARG" ;;
     h) echo ""
         echo "Version 2.0: released on June 10th, 2013"
         echo ""
         echo "Usage: ./ASDEG.sh -p <PATH> [options ... ]"
         echo "  -p <PATH>  Full path of data directory <EX. ~/SNP0.05_500/data_1M_5_5_0_20>"
         echo "  -S <PATH>  Full path of script directory < EX.~/SNP0.05_500/data_1M_5_5_0_20>"
         echo "  -G <file_name> snp2gene file . "
         echo "  -t <int> Number of thread [Default 1]"
         echo "  -c <file_name> Case file name. If more than two are provided, use case1,case2,...,casen"
         echo "  -C <file_name> Control file name. If more than two are provided, use control1,control2,...,controln"
         echo "  -g <int> the column where the gene name locate"
         echo "  -o <DIR> the output directory [Default ./ASDEG]"
         echo "  -h         help"
         echo ""
         exit;
         ;;
     esac
 done

###建output資料夾

output=./$output_dir

if ! [ -d $output ];then

  mkdir $output

fi

##############################Use for Creating tabix index #########################################

echo "[executing program] create index "

echo "${case}" > $output/case_tmp
case_file=($(less $output/case_tmp | sed 's/,/\n/g'))
case_count=`less $output/case_tmp | sed 's/,/\n/g' | wc -l`

echo "${control}" > $output/control_tmp
control_file=($(less $output/control_tmp | sed 's/,/\n/g'))
control_count=`less $output/control_tmp | sed 's/,/\n/g' | wc -l`

for (( i=0; i<$case_count; i++))

    do
       if [ -e ${data_path}/${case_file[i]}.gz -a ${data_path}/${case_file[i]}.gz.tbi ];then

            continue
       else
            ${script_path}/create_index.sh ${data_path}/${case_file[i]} $script_path
    fi
done

for (( i=0; i<$control_count; i++))

    do
       if [ -e ${data_path}/${control_file[i]}.gz -a ${data_path}/${control_file[i]}.gz.tbi ];then

            continue
       else
            ${script_path}/create_index.sh ${data_path}/${control_file[i]} $script_path
    fi
done

date '+%T'

rm $output/case_tmp $output/control_tmp

############################project開始################################################

echo "[executing program] use gene name to split snp2gene file "

date '+%T'

$script_path/split_gene_name $data_path $snp2gene $gene_name_column $script_path $output

trimming_path=$output/${snp2gene}_new_folder

a=`less $trimming_path/${snp2gene}_new | grep -v 'pos' | wc -l`
a=`expr $a \* 2 `

echo "[executing program] data processing (trimming fasta)"

$script_path/Trimming_multi.sh $data_path $trimming_path $thread $read_length $case $control $script_path $output

b=`less $trimming_path/log | wc -l`

while [ $a -gt $b ]

   do

      b=`less $trimming_path/log | wc -l`

      continue
done

date '+%T'

echo "[executing program] data processing finished"

echo "[executing program] Fisher Exact Test and Combine P_value"

date '+%T'

$script_path/Fisher_Test.sh $data_path $trimming_path ${snp2gene}_new $case $control $script_path $output $percentage

date '+%T'

echo "[executing program] FDR Control by BH"

$script_path/FDR_Control.sh $output $script_path

date '+%T'

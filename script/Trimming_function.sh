#!/bin/bash

##### Setting parameters #####
data_path=$1
trimming_path=$2
script_path=$3
read_length=$4
sample=$5
snp_chromsome=$6
snp_gene=$7
snp_coordinate_max=$8
gene_file=$9
output=${10}
################################
echo "${sample}" >  $trimming_path/"${snp_gene}"_folder/sample_tmp_$8_$5
sample_file=($(less $trimming_path/"${snp_gene}"_folder/sample_tmp_$8_$5 | sed 's/,/\n/g' | cut -d '.' -f 1))
sample_count=`less $trimming_path/"${snp_gene}"_folder/sample_tmp_$8_$5 | sed 's/,/\n/g' | wc -l`

    for (( i=0; i<$sample_count; i=i+1))

         do
            k=`expr ${snp_coordinate_max} ` ##k is snp position
            
            ###xMxS && xM
            
	    $script_path/tabix ${data_path}/"${sample_file[i]}".sam.gz "${snp_chromsome}":$k-$k | grep -v -G "[0-9]*M[0-9]*N[0-9]*M" | grep -v -G "[0-9]*S[0-9]*M" | cut -f 2,4,6,10 | grep -v 'H' | grep -v 'I' | grep -v 'D'| sed 's/M/\t/g' | sed 's/S/\tS/g' | awk '{if('$k'-$2+1<$3 && '$k'-$2+1>0){$1='$k'-$2+1}else{$1="*"}}{print$1"\t"$2"\t"$4"\t"$5"\t"$6}' | awk '{if($4!="S" && $1!="*"){print$1"\t"$3}else if($4=="S" && $1!="*"){print$1"\t"$5}}' > $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test

            ### xSxM && xSxMxS
    
            $script_path/tabix ${data_path}/"${sample_file[i]}".sam.gz "${snp_chromsome}":$k-$k | grep -G "[0-9]*S[0-9]*M" | grep -v -G "[0-9]*S[0-9]*M[0-9]*N" | cut -f 2,4,6,10 | grep -v 'H' | grep -v 'I' | grep -v 'D'| sed 's/M/\t/g' | sed 's/S/\tS\t/g' | awk '{if('$k'-$2-$3+1<$5 && '$k'-$2-$3+1>0){$1='$k'-$2-$3+1}else{$1="*"}}{print$1"\t"$2"\t"$6"\t"$7"\t"$8}' | awk '{if($4!="S" && $1!="*"){print$1"\t"$3}else if($4=="S" && $1!="*"){print$1"\t"$5}}' >> $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test

            ### xMxNxM && xMxNxMxS

            $script_path/tabix ${data_path}/"${sample_file[i]}".sam.gz "${snp_chromsome}":$k-$k | grep -G "[0-9]*M[0-9]*N[0-9]*M" | grep -v -G "[0-9]*S[0-9]*M[0-9]*N[0-9]*M" | grep -v -G "[0-9]*M[0-9]*N[0-9]*M[0-9]*N[0-9]*M" | cut -f 2,4,6,10 | grep -v 'H' | grep -v 'I' | grep -v 'D'| sed 's/M/\t/g' | sed 's/N/\t/g' | sed 's/S/\tS/g' | awk '{if('$k'-$2+1>0 && '$k'-$2+1<=$3){$1='$k'-$2+1}else if('$k'-$2-$3-$4+1>0 && '$k'-$2-$3-$4+1<=$5){$1='$k'-$2-$4+1}else{$1="*"}}{print$1"\t"$2"\t"$6"\t"$7"\t"$8}'  | awk '{if($4!="S" && $1!="*"){print$1"\t"$3}else if($4=="S" && $1!="*"){print$1"\t"$5}}' >> $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test

            ### xSxMxNxM && xSxMxNxMxS

	    $script_path/tabix ${data_path}/"${sample_file[i]}".sam.gz "${snp_chromsome}":$k-$k | grep -G "[0-9]*S[0-9]*M[0-9]*N[0-9]*M" | grep -v -G "[0-9]*S[0-9]*M[0-9]*N[0-9]*M[0-9]*N" | cut -f 2,4,6,10 | grep -v 'H' | grep -v 'I' | grep -v 'D'| sed 's/M/\t/g' | sed 's/N/\t/g' | sed 's/S/\tS\t/g' | awk '{if('$k'-$2-$3+1>0 && '$k'-$2-$3+1<=$5){$1='$k'-$2-$3+1}else if('$k'-$2-$3-$5-$6+1>0 && '$k'-$2-$3-$5-$6+1<=$7){$1='$k'-$2-$3-$6+1}else{$1="*"}}{print$1"\t"$2"\t"$8"\t"$9"\t"$10}' | awk '{if($4!="S" && $1!="*"){print$1"\t"$3}else if($4=="S" && $1!="*"){print$1"\t"$5}}' >> $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test

            read_count=`less $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test | wc -l`

            if [ "$read_count" -ne 0 ]; then

            trimming_sample=($(less $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test | cut -f 1 | sort | uniq ))
            trimming_count=`less $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test | cut -f 1 | sort | uniq | wc -l`

            for (( m=0; m<$trimming_count; m=m+1 ))

              do
                
                length=`expr ${trimming_sample[m]} `
              
                less $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test | awk '$1=='$length'' | cut -f 2 | cut -c $length >> $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}.trimming

          done

      fi
               echo "${snp_gene} finished " >> $trimming_path/log
               rm $trimming_path/"${snp_gene}"_folder/"${sample_file[i]}"_"${snp_gene}"_${k}_test
                        
  done

   rm $trimming_path/"${snp_gene}"_folder/sample_tmp_$8_$5

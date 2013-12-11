#!/bin/bash

########## Setting parameters ##########
main_path=$1
trimming_path=$2
gene_snp_file=$3
case=$4
control=$5
script_path=$6
output=$7
percentage=$8
########################################
##########create folder#################
if [ -d $trimming_path/count_result_s ];then

   rm -r $trimming_path/count_result_s
   mkdir $trimming_path/count_result_s

else

   mkdir $trimming_path/count_result_s
fi

final_result=$trimming_path/count_result_s

########################################
gene_name=($(ls $trimming_path | grep "gene_name" | cut -d '_' -f 3 ))
total_gene=`ls $trimming_path | grep "gene_name" | cut -d '_' -f 3 | wc -l`

echo "${case}" > $final_result/case_tmp
case_file=($(less $final_result/case_tmp | sed 's/,/\n/g' | cut -d '.' -f 1))
case_count=`less $final_result/case_tmp | sed 's/,/\n/g' | wc -l`

echo "${control}" > $final_result/control_tmp
control_file=($(less $final_result/control_tmp | sed 's/,/\n/g' | cut -d '.' -f 1))
control_count=`less $final_result/control_tmp | sed 's/,/\n/g' | wc -l`

total_sample=`expr $case_count + $control_count `

for ((i=0; i<$total_gene; i=i+1))   #####總共有多少個基因需要處理

  do
     gene_tag=0
     gene_chromosome=($(less $trimming_path/gene_name_${gene_name[i]} | cut -f 1 | uniq))
     fix=`ls $trimming_path/${gene_name[i]}_folder | grep 'trimming' | wc -l`  ###判斷資料夾裡面是否有資料需要處理
 
     if [ $fix -eq 0 ];then   ###若是沒有，則不處理

        continue
     else
        
     snp_position=($(less $trimming_path/gene_name_${gene_name[i]} | cut -f 2 | grep -v 'pos' | sort | uniq))
     snp_count=`less $trimming_path/gene_name_${gene_name[i]} | cut -f 2 | grep -v 'pos' | sort | uniq | wc -l`

    for (( j=0; j<$snp_count; j++ ))

     do
        snp_tag=0
        sample_count=`ls $trimming_path/${gene_name[i]}_folder | grep 'trimming' | grep "_${snp_position[j]}." | wc -l`

        if [ $sample_count -ne $total_sample ];then ###若是需要處理所有sample和原本的sample數不相同，則不處理

           continue

        else
            
           allele_case=($(cat $trimming_path/${gene_name[i]}_folder/*_${snp_position[j]}.trimming | sort | uniq))
           allele_count=`cat $trimming_path/${gene_name[i]}_folder/*_${snp_position[j]}.trimming | sort | uniq | wc -l`

#####此部分的程式是為了判斷case和control表現量最高的allele是否是同一個，若是兩者表現較高的allele是相同的，則不做後續統計測試
         
                for (( case=0; case<$case_count; case++ ))
                   do
                      less $trimming_path/${gene_name[i]}_folder/${case_file[case]}_${gene_name[i]}_${snp_position[j]}.trimming | sort | uniq -c | sort -nrk 1 | head -n 1 | awk '{print$2}' >> $final_result/${gene_name[i]}_${snp_position[j]}_case_obvious_allele

                done

                for (( control=0; control<$control_count; control++ ))
                   do
                      less $trimming_path/${gene_name[i]}_folder/${control_file[control]}_${gene_name[i]}_${snp_position[j]}.trimming | sort | uniq -c | sort -nrk 1 | head -n 1 | awk '{print$2}' >> $final_result/${gene_name[i]}_${snp_position[j]}_control_obvious_allele

                   done
                
                treatment_allele_case=`less $final_result/${gene_name[i]}_${snp_position[j]}_case_obvious_allele | sort | uniq | wc -l`
                control_allele_case=`less $final_result/${gene_name[i]}_${snp_position[j]}_control_obvious_allele | sort | uniq | wc -l` 

                if [ $treatment_allele_case -eq 1 -a $control_allele_case -eq 1 ];then
                   
                   merge_allele_case=`cat $final_result/${gene_name[i]}_${snp_position[j]}_*_obvious_allele | sort | uniq | wc -l`
                   
                   if [ $merge_allele_case -eq 1 ];then
 
                       echo "${snp_position[j]}" >> $final_result/nodo

                   else ###這一部分為需要做Fisher Exact Test的部分

                      if [ $gene_tag -eq 0 ];then
                         echo "=============================" >> $final_result/${gene_name[i]}_information
                         echo "Gene_Chromosome:${gene_chromosome[0]}" >> $final_result/${gene_name[i]}_information
                         echo "=============================" >> $final_result/${gene_name[i]}_information
                         gene_tag=1
                      fi

                      if [ $snp_tag -eq 0 ];then
                         echo "***SNP Position:${snp_position[j]}***" >> $final_result/${gene_name[i]}_information
                         snp_tag=1
                      fi          
                      
                      for (( case=0; case<$case_count; case++ ))

                        do

                          for (( control=0; control<$control_count; control++ ))
                            
                            do
                               echo "alleles <- matrix(c(" >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R

                               case_total_count=0
                               control_total_count=0

                               for (( al=0; al<$allele_count; al++ ))

                                 do
                                    case_allele_read_count=`less $trimming_path/${gene_name[i]}_folder/${case_file[case]}_${gene_name[i]}_${snp_position[j]}.trimming | grep -w "${allele_case[al]}" | wc -l`
                                    control_allele_read_count=`less $trimming_path/${gene_name[i]}_folder/${control_file[control]}_${gene_name[i]}_${snp_position[j]}.trimming | grep -w "${allele_case[al]}" | wc -l`
                                    
                                    echo -e "allele:${allele_case[al]}\tread_number:$case_allele_read_count" >> $final_result/${gene_name[i]}_${snp_position[j]}_case
                                    echo -e "allele:${allele_case[al]}\tread_number:$control_allele_read_count" >> $final_result/${gene_name[i]}_${snp_position[j]}_control

                                    if [ $al -le `expr $allele_count - 2` ]; then    ###輸入呼叫R做Fisher Exact Test的格式並做Fisher Exact Test

                                       echo "$case_allele_read_count,$control_allele_read_count," >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R
                                    else
                                       echo "$case_allele_read_count,$control_allele_read_count),nr=2)" >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R
                                       echo -e "sink(\"$final_result/${gene_name[i]}_${snp_position[j]}_$case$control.txt\")" >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R
                                       echo "fisher.test(alleles)" >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R
                                       echo "sink()" >> $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R
                                       $script_path/R --no-save < $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R >> $final_result/R.log
                                       rm $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.R                                       
                                    fi    ###做完Fisher Exact Test
                                 done   ###allele count 的 loop
                                 echo "Case Sample=${case_file[case]}" >> $final_result/${gene_name[i]}_information
                                 cat $final_result/${gene_name[i]}_${snp_position[j]}_case >> $final_result/${gene_name[i]}_information
                                 echo "Control Sample=${control_file[control]}" >> $final_result/${gene_name[i]}_information
                                 cat $final_result/${gene_name[i]}_${snp_position[j]}_control >> $final_result/${gene_name[i]}_information
                                 p_value=`less $final_result/${gene_name[i]}_${snp_position[j]}_$case$control.txt | grep 'p-value' | sed 's/</_/g' | sed 's/=/_/g' | cut -d '_' -f 2`
                                 echo "${p_value}" >> $final_result/${gene_name[i]}_FDR
                                 echo -e "[p_value]\t${p_value}" >> $final_result/${gene_name[i]}_information
                                 echo "------------------------------" >> $final_result/${gene_name[i]}_information
                                 rm $final_result/${gene_name[i]}_${snp_position[j]}_case $final_result/${gene_name[i]}_${snp_position[j]}_control          
                             done   ###control file 的loop
                        done   ###case file 的 loop                          
                      fi   ###case和control group表現量較高的allele是不相同的才做fisher exact test
               else

                   #rm $trimming_path/count_result_s/${gene_name[i]}_${snp_position[j]}_*_obvious_allele
                    echo "nodo" > tmp_nodo
                    rm tmp_nodo

               fi ###case group and control group各自表現量較高的allele是相同的才做處理

       fi   ###若是需要處理所有sample和原本的sample數不相同，則不處理
          
   done   ### snp loop

   rm $final_result/${gene_name[i]}_*_obvious_allele

###FDR正式開始###

if [ -f $final_result/${gene_name[i]}_FDR ];then

   FDR_pvalue=($(less $final_result/${gene_name[i]}_FDR | cut -f 1 ))
   FDR_pvalue_count=`less $final_result/${gene_name[i]}_FDR | cut -f 1 | wc -l`
   
   for (( f=0; f<$FDR_pvalue_count; f++ ))

       do
          if [ $f -eq 0 ];then

              echo "x <- c(" >> $final_result/${gene_name[i]}_FDR.R
          fi
          if ! [ $f -eq `expr $FDR_pvalue_count - 1` ]; then
        
             echo "${FDR_pvalue[f]}," >> $final_result/${gene_name[i]}_FDR.R

          else
           
             echo "${FDR_pvalue[f]})" >> $final_result/${gene_name[i]}_FDR.R
             echo -e "sink(\"$final_result/${gene_name[i]}_SNP_FDR\")" >> $final_result/${gene_name[i]}_FDR.R
             echo -e "p.adjust(x, method = \"BH\", n = length(x))" >> $final_result/${gene_name[i]}_FDR.R
             echo "sink()" >> $final_result/${gene_name[i]}_FDR.R
             $script_path/R --no-save < $final_result/${gene_name[i]}_FDR.R >> $final_result/R.log
             rm -r $final_result/${gene_name[i]}_FDR.R
             less $final_result/${gene_name[i]}_SNP_FDR | sed 's/] /]/g' | cut -d ']' -f 2 | sed 's/ /\n/g' > $final_result/tmp
             mv $final_result/tmp $final_result/${gene_name[i]}_SNP_FDR
             rm $final_result/R.log       
          fi
       done
    combine_p_value=($( less $final_result/${gene_name[i]}_SNP_FDR | cut -f 1 ))
    combine_p_value_count=`less $final_result/${gene_name[i]}_SNP_FDR | cut -f 1 | wc -l`
  
    for (( p=0; p<$combine_p_value_count; p++ ))

        do
           $script_path/correct_pvalue ${combine_p_value[p]} > $final_result/correct
           correct_pvalue=`less $final_result/correct | cut -d '=' -f 2`

           if [ $p -eq 0 ];then

              echo "x <- c(" >> $final_result/${gene_name[i]}_combine.R
              IUT_p_value=`expr $correct_pvalue`
              correct_rev=$correct_pvalue
           
           else
            
              if [ $(echo "$IUT_p_value < $correct_pvalue " | bc) -eq 1 ];then
 
              IUT_p_value=`expr $correct_pvalue `
              correct_rev=$correct_pvalue       
              
              fi
           fi
           

           if ! [ $p -eq `expr $combine_p_value_count - 1` ]; then
 
              echo "${combine_p_value[p]}," >> $final_result/${gene_name[i]}_combine.R
    
           else
          
              echo "${combine_p_value[p]})" >> $final_result/${gene_name[i]}_combine.R
              echo -e "sink(\"$final_result/${gene_name[i]}_combine\")" >> $final_result/${gene_name[i]}_combine.R
              echo "pnorm(sum(qnorm(x))/sqrt(length(x)))" >> $final_result/"${gene_name[i]}"_combine.R
              echo "pchisq(sum(-2*log(x)),2*length(x),low=F)" >> $final_result/"${gene_name[i]}"_combine.R
              echo "sink()" >> $final_result/"${gene_name[i]}"_combine.R
              $script_path/R --no-save < $final_result/"${gene_name[i]}"_combine.R >> $final_result/R.log
              rm $final_result/"${gene_name[i]}"_combine.R
              combine_pvalue=($(less $final_result/${gene_name[i]}_combine | cut -d ']' -f 2))
              echo -e "gene_name:${gene_name[i]}\t${combine_pvalue[0]}" >> $output/stoufferst_result
              echo -e "gene_name:${gene_name[i]}\t${combine_pvalue[1]}" >> $output/fisherst_result      
           fi
        done
       $script_path/new_pvalue $correct_rev > $final_result/new
       IUT_p_value=`less $final_result/new`        
       echo -e "gene_name:${gene_name[i]}\t$IUT_p_value" >> $output/IUT_result
       rm $final_result/new $final_result/correct $final_result/${gene_name[i]}_*.txt $final_result/${gene_name[i]}_SNP_FDR $final_result/${gene_name[i]}_combine
fi

 fi  ###判斷資料夾裡面是否有資料需要處理
done   ###基因的迴圈

rm $final_result/case_tmp $final_result/control_tmp $final_result/nodo 

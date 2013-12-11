#!/bin/bash

##### Setting parameters #####
data_path=$1
script_path=$2
##############################

file_name=($(ls $data_path | grep '_result'))
method=($(ls $data_path | grep '_result' | cut -d '_' -f 1))
file_count=`ls $data_path | grep '_result' | wc -l`

for (( f=0; f<$file_count; f=f+1 ))

do
   gene_name=($(less $data_path/${file_name[f]} | awk '$3=1-$2' | sed 's/ /\t/g' | awk '$2 < 1' | sort -n -k 2 | cut -f 1 | cut -d ':' -f 2))
   pvalue=($(less $data_path/${file_name[f]} |  awk '$3=1-$2' | sed 's/ /\t/g' | awk '$2 < 1' | sort -n -k 2 | cut -f 2 ))
   gene_count=`less $data_path/${file_name[f]} |  awk '$3=1-$2' | sed 's/ /\t/g' | awk '$2 < 1' | sort -n -k 2 | cut -f 2 | wc -l`
   
#   echo $gene_count
   for (( j=0; j<$gene_count; j=j+1 ))

    do
       if [ $j -eq 0 ];then

         echo "pvalue <- c(" >> $data_path/${method[f]}.R
         echo "pgene <- c(" >> $data_path/${method[f]}.gene
       fi

       if [ $j -le `expr $gene_count - 2` ]; then

         echo "${pvalue[j]}," >> $data_path/${method[f]}.R
         echo "(${gene_name[j]})," > $data_path/tmp
         less $data_path/tmp | sed 's/(/"/g' | sed 's/)/"/g' >> $data_path/${method[f]}.gene
         rm $data_path/tmp        

       else
         
        echo "${pvalue[j]})" >> $data_path/${method[f]}.R
        echo "(${gene_name[j]})" > $data_path/tmp
        less $data_path/tmp | sed 's/(/"/g' | sed 's/)/")/g' >> $data_path/${method[f]}.gene
        cat $data_path/${method[f]}.gene >> $data_path/${method[f]}.R
        rm $data_path/tmp $data_path/${method[f]}.gene
        echo "sorted.pvalue <- sort(pvalue)" >> $data_path/${method[f]}.R
        echo "j.alpha <- (1:length(pvalue))*(.05/length(pvalue))" >> $data_path/${method[f]}.R
        echo "dif <- sorted.pvalue-j.alpha" >> $data_path/${method[f]}.R
        echo "neg.dif<-dif[dif<0]" >> $data_path/${method[f]}.R
        echo "pos.dif<-neg.dif[length(neg.dif)]">> $data_path/${method[f]}.R
        echo "index<-dif==pos.dif" >> $data_path/${method[f]}.R
        echo "p.cutoff<-sorted.pvalue[index]" >> $data_path/${method[f]}.R
        echo "p.sig<-pvalue[pvalue<=p.cutoff]" >> $data_path/${method[f]}.R
        echo "p.gene<-pgene[pvalue<=p.cutoff]" >> $data_path/${method[f]}.R
        echo "sink($data_path/FDR_correct)" > $data_path/temp
        less $data_path/temp | sed 's/sink(/sink("/g' | sed 's/correct)/correct")/g' >> $data_path/${method[f]}.R
        echo "p.sig" >> $data_path/${method[f]}.R
        echo "p.gene" >> $data_path/${method[f]}.R
        echo "sink()" >> $data_path/${method[f]}.R
        $script_path/R --no-save < $data_path/${method[f]}.R >> $data_path/R.log
        rm $data_path/${method[f]}.R $data_path/temp
       fi
   done
    
    reject_pvalue=($(less $data_path/FDR_correct | grep -v "numeric" | grep -v "character" | grep -v '"' |  cut -d ']' -f 2 | colrm 1 1 | sed 's/  /\n/g' |sed 's/ /\n/g'))
    reject_gene=($(less $data_path/FDR_correct | grep -v "character" | grep -v "numeric" | grep '"' | sed 's/"//g' | cut -d ']' -f 2 | colrm 1 1 | sed 's/  /\n/g' |sed 's/ /\n/g'))
    reject_count=`less $data_path/FDR_correct | grep -v "numeric" | grep -v "character" |  grep -v '"' |  cut -d ']' -f 2 | colrm 1 1 | sed 's/  /\n/g' |sed 's/ /\n/g' | wc -l`

#    echo "$reject_count"

    if [ $reject_count -gt 0 ];then

      echo "##################################" >> $data_path/${method[f]}_FDR_final
      echo "method:${method[f]}" >> $data_path/${method[f]}_FDR_final
      echo "reject HO gene count:$reject_count" >> $data_path/${method[f]}_FDR_final
      echo "##################################" >> $data_path/${method[f]}_FDR_final
      for (( re=0; re<$reject_count; re=re+1 ))
      do
         echo -e "${reject_gene[re]}\t${reject_pvalue[re]}" >> $data_path/${method[f]}_FDR_final
   
     done
    else
        echo "##################################" >> $data_path/${method[f]}_FDR_final
        echo "method:${method[f]}" >> $data_path/${method[f]}_FDR_final
        echo "##################################" >> $data_path/${method[f]}_FDR_final
        echo "Do not have any reject H0 gene " >> $data_path/${method[f]}_FDR_final

    fi

    rm $data_path/FDR_correct
        
done

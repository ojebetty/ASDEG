ASDEG
=====

ASDEG : Allele Specific Differentail Expression Gene

The following is ASDEG usage. 

Version 2.0: released on June 10th, 2013

Usage: ./ASDEG.sh -p < PATH > [options ... ]

  -p < PATH >  Full path of data directory < EX. ~/SNP0.05_500/data_1M_5_5_0_20 

  -S < PATH >  Full path of script directory < EX.~/SNP0.05_500/data_1M_5_5_0_20>

  -G < file_name > snp2gene file .

  -t < int > Number of thread [Default 1]

  -c < file_name > Case file name. If more than two are provided, use case1,case2,...,casen

  -C < file_name > Control file name. If more than two are provided, use control1,control2,...,controln

  -g < int > the column where the gene name locate

  -o < DIR > the output directory [Default ./ASDEG]

  -h         help

#!/bin/bash

##########資料處理################


#less $1 | sed 's/-/_/g' > tmp
#mv tmp $1

######### use tabix to create index ##

$2/bgzip $1
$2/tabix -p sam $1.gz


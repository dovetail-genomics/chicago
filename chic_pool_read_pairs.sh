#!/bin/bash

if [ $# -lt 1 ]
   then
      echo "Usage: chic_pool_read_pairs.sh <chic_split_by_chr output file>"
      exit
fi

echo "baitID	otherEndID	N	otherEndLen	distSign"; 

awk '{ 
    if (!baitOtherEndN[$1"\t"$2]){ 
         baitOtherEndN[$1"\t"$2] = 1; 
         baitOtherEndInfo[$1"\t"$2] = $3"\t"$4 
    }
    else{ 
         baitOtherEndN[$1"\t"$2]++; 
    }  
}END{ 
    for (key in baitOtherEndN){
         print key"\t"baitOtherEndN[key]"\t"baitOtherEndInfo[key];
    }
}' $1 | sort -k1,1 -k2,2n

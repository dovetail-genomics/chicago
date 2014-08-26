#!/bin/bash

if [ $# -lt 4 ]
  then
    echo "Usage: process_chic_single_core.sh <bamfile> <baits-with-fragment-ids-file> <digest-with-fragment-ids-file> <sample-name>"
    exit
fi

pipelinedir="/bi/group/sysgen/mikhail/my_chic/pipeline"

processhicup="${pipelinedir}/process_hicup.sh"
chicnosplit="${pipelinedir}/chic_no_split.sh"
chicpoolreadpairs="${pipelinedir}/chic_pool_read_pairs.sh"

bam=$1
bam0=`basename ${bam}`
bamname=`echo ${bam0%\.*}`
baitfendsid=$2
digestbed=$3
samplename=$4

mkdir -p sample_${samplename}

echo "Processing sample ${samplename}..."
echo "Using bam file ${bam}"
echo "Using bait map file ${baitfendsid}"
echo "Using fragment map file ${digestbed}"

${processhicup} ${bam} ${baitfendsid} ${digestbed} sample_${samplename}

echo "Pooling read pairs..."

echo "baitID	otherEndID	N	otherEndLen	distSign" > sample_${samplename}/${samplename}_bait_otherEnd_N_len_distSign.txt
awk '{ 
    if (!baitOtherEndN[$14"\t"$18]){ 
         baitOtherEndN[$14"\t"$18] = 1; 
         baitOtherEndInfo[$14"\t"$18] = $20"\t"$21 
    }
    else{ 
         baitOtherEndN[$14"\t"$18]++; 
    }  
}END{ 
    for (key in baitOtherEndN){
         print key"\t"baitOtherEndN[key]"\t"baitOtherEndInfo[key];
    }
}' sample_${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe | sort -k1,1 -k2,2n >> sample_${samplename}/${samplename}_bait_otherEnd_N_len_distSign.txt


echo "Done! The final output file with pooled reads per interaction is sample_${samplename}/${samplename}_bait_otherEnd_N_len_distSign.txt"

#!/bin/bash

if [ $# -lt 4 ]
  then
    echo "Usage: bam2chicago.sh <bamfile> <baitmap-file> <digest-rmap-file> <sample-name> [nodelete]"
    exit
fi

SCRIPT=$(readlink -f $0)
pipelinedir=`dirname $SCRIPT`

processhicup="${pipelinedir}/process_hicup.sh"
chicnosplit="${pipelinedir}/chic_no_split.sh"
chicpoolreadpairs="${pipelinedir}/chic_pool_read_pairs.sh"

bam=$1
bam0=`basename ${bam}`
bamname=`echo ${bam0%\.*}`
baitfendsid=$2
digestbed=$3
samplename=$4
nodelete=$5

mkdir -p ${samplename}

echo "Processing sample ${samplename}..."
echo "Using bam file ${bam}"
echo "Using bait map file ${baitfendsid}"
echo "Using fragment map file ${digestbed}"

${processhicup} ${bam} ${baitfendsid} ${digestbed} ${samplename}

echo "Pooling read pairs..."

echo "baitID	otherEndID	N	otherEndLen	distSign" > ${samplename}/${samplename}_bait_otherEnd_N_len_distSign.txt
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
}' ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe | sort -k1,1 -k2,2n -T ${samplename} >> ${samplename}/${samplename}.chinput

mv ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fless06_withDistSignLen.bedpe ${samplename}/${bamname}_ambiguous_alignments.bedpe

if [ "$nodelete" != "nodelete" ]; then
	echo "Removing intermediate files..."
	rm ${samplename}/${bamname}_mappedToBaits_baitOnRight.bedpe
	rm ${samplename}/${bamname}_mappedToBaits.bedpe
	rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe
	rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe
	rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe
fi

echo "Done! The file to be used for Chicago R package input is ${samplename}/${samplename}.chinput"

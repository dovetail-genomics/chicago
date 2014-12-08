#!/bin/bash

if [ $# -lt 1 ]
   then
      echo "Usage: chic_no_split.sh <hicup_process.sh output file> <sample-id>"
      exit
fi

inputfile=$1
sampleid=$2

awk -v sampleid=$sampleid 'BEGIN{OFS="\t"}{ print $14,$18,$20,$21 }' $inputfile > ${sampleid}_bait_otherEnd_len_distSign.txt


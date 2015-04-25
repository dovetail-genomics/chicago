#!/bin/bash

if [ $# -lt 3 ]
  then
    echo "Usage: process_hicup.sh <bamfile> <baitmap-file> <rmap-file> [<output-dir>]"
    exit
fi

bam=$1
bamname=`echo ${bam%\.*}`

baitfendsid=$2
digestbed=$3

outputdir=.

if [ $# -gt 3 ]
  then
    outputdir=$4
fi

#echo "Converting the input bam file ${bam} to a paired-end bed file ${bamname}.bedpe..."
#bedtools bamtobed -bedpe -i $bam > ${outputdir}/${bamname}.bedpe

baitcolno=`awk '{print NF; exit}' ${baitfendsid}`

if [ $baitcolno -gt 4 ]
  then
	echo "Baitfendsid file contains >4 columns. Probably means that full bait restriction fragment map was provided instead of valid fends, but it's ok for HiCUP output at least. Trimming the extra columns and saving the result in ${baitfendsid}_4col.txt..."
    if ! [ -e ${baitfendsid}_4col.txt ] 
		then 
			awk '{print $1"\t"$2"\t"$3"\t"$4}' ${baitfendsid} > ${baitfendsid}_4col.txt
	fi
    baitfendsid=`echo ${baitfendsid}_4col.txt`
fi

echo "Intersecting with bait fragments (using min overhang of 0.6)..."
#bedtools pairtobed -a ${outputdir}/${bamname}.bedpe -b $baitfendsid -f 0.6 > ${outputdir}/${bamname}_mappedToBaits.bedpe
bedtools pairtobed -abam $bam -bedpe -b $baitfendsid -f 0.6 > ${outputdir}/${bamname}_mappedToBaits.bedpe

echo "Flipping all reads that overlap with the bait on to the right-hand side..."
#awk 'BEGIN{OFS="\t"}{if($1==$11 && $3>=$12 && $2<=$13){print $4,$5,$6,$1,$2,$3,$7,$8,$10,$9,$11,$12,$13,$14,$15}else{print $0}}' ${outputdir}/${bamname}_mappedToBaits.bedpe > ${outputdir}/${bamname}_mappedToBaits_baitOnRight.bedpe
awk 'BEGIN{OFS="\t"}{minRight=$13<$3?$13:$3; maxLeft=$12>$2?$12:$2; if($1==$11 && (minRight-maxLeft)/($3-$2)>=0.6){print $4,$5,$6,$1,$2,$3,$7,$8,$10,$9,$11,$12,$13,$14,$15}else{print $0}}' ${outputdir}/${bamname}_mappedToBaits.bedpe > ${outputdir}/${bamname}_mappedToBaits_baitOnRight.bedpe

echo "Intersecting with bait fragments again to produce a list of bait-to-bait interactions that can be used separately; note they will also be retained in the main output..."
bedtools intersect -a ${outputdir}/${bamname}_mappedToBaits_baitOnRight.bedpe -wo -f 0.6 -b $baitfendsid > ${outputdir}/${bamname}_bait2bait.bedpe

#echo "Removing reads that failed the min overhang filter..."
#awk '{if ($0!~/\-1\t\-1/){print $0}}' ${outputdir}/${bamname}_bait2bait.bedpe > ${outputdir}/${bamname}_bait2bait_fmore06.bedpe

echo "Intersecting with restriction fragments (using min overhang of 0.6)..."
bedtools intersect -a ${outputdir}/${bamname}_mappedToBaits_baitOnRight.bedpe -wao -f 0.6 -b $digestbed > ${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe

echo "Removing reads that failed the min overhang filter (saving separately into ${bamname}_mappedToBaitsBoRAndRFrag_fless06.bedpe)..."
awk -v fless=${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag_fless06.bedpe '{if ($0~/\-1\t\-1/){print $0 > fless} else {print $0}}' ${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe > ${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe

echo "Adding frag length and signed distance from bait; removing self-ligation fragments (if any; not expected with HiCUP input)..."

perl -ne '{ chomp $_; my @a = split /\t/, $_; my $l = $a[16]-$a[15]; my $d = "NA"; if ($a[10] eq $a[14]){ my $midB = $a[11]+int(($a[12]-$a[11])/2+0.5); my $midR = $a[15]+int(($a[16]-$a[15])/2+0.5); $d = $midR - $midB; } if ($d ne "0") { print "$_\t$l\t$d\n" }}' ${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe > ${outputdir}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe

echo "The file with reads aligned to baits and other ends is in the ${outputdir} folder and is called ${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe"

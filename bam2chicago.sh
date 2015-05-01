#!/bin/bash

if [ $# -lt 4 ]
  then
    echo "Usage: bam2chicago.sh <bamfile> <baitmap-file> <digest-rmap-file> <sample-name> [nodelete]"
    exit
fi

command -v bedtools >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute bedtools. Check that it's installed and added to PATH. Aborting."; exit 1; }
command -v awk >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute awk. Aboring."; exit 1; }
command -v perl >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute perl. Aboring."; exit 1; }

#SCRIPT=$(readlink -f $0)
#pipelinedir=`dirname $SCRIPT`

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
echo "Using baitmap file ${baitfendsid}"
echo "Using digest map (rmap) file ${digestbed}"

baitcolno=`awk '{print NF; exit}' ${baitfendsid}`

if [ $baitcolno -gt 4 ]; then
    echo "Baitmap file contains >4 columns. Checking if ${baitfendsid}_4col.txt exists..."
    if ! [ -e ${baitfendsid}_4col.txt ]; then
	   echo "It doesn't. So trimming the extra columns and saving the result in ${baitfendsid}_4col.txt..."
           awk '{print $1"\t"$2"\t"$3"\t"$4}' ${baitfendsid} > ${baitfendsid}_4col.txt
    else
           echo "Found ${baitfendsid}_4col.txt"
    fi
    baitfendsid=`echo ${baitfendsid}_4col.txt`
fi

echo "Intersecting with bait fragments (using min overhang of 0.6)..."
bedtools pairtobed -abam $bam -bedpe -b $baitfendsid -f 0.6 > ${samplename}/${bamname}_mappedToBaits.bedpe

echo "Flipping all reads that overlap with the bait on to the right-hand side..."
awk 'BEGIN{ OFS="\t" }
     {
     minRight=$13<$3?$13:$3; 
     maxLeft=$12>$2?$12:$2; 
     if($1==$11 && (minRight-maxLeft)/($3-$2)>=0.6){
              print $4,$5,$6,$1,$2,$3,$7,$8,$10,$9,$11,$12,$13,$14,$15
     }
     else {
              print $0
     }
}' ${samplename}/${bamname}_mappedToBaits.bedpe > ${samplename}/${bamname}_mappedToBaits_baitOnRight.bedpe

if [ "$nodelete" != "nodelete" ]; then
	rm ${samplename}/${bamname}_mappedToBaits.bedpe
fi

echo "Intersecting with bait fragments again to produce a list of bait-to-bait interactions that can be used separately; note they will also be retained in the main output..."
bedtools intersect -a ${samplename}/${bamname}_mappedToBaits_baitOnRight.bedpe -wo -f 0.6 -b $baitfendsid > ${samplename}/${bamname}_bait2bait.bedpe

echo "Intersecting with restriction fragments (using min overhang of 0.6)..."
bedtools intersect -a ${samplename}/${bamname}_mappedToBaits_baitOnRight.bedpe -wao -f 0.6 -b $digestbed > ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe

if [ "$nodelete" != "nodelete" ]; then
	rm ${samplename}/${bamname}_mappedToBaits_baitOnRight.bedpe
        lessfilename=/dev/null
else
        lessfilename=${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fless06.bedpe
fi

echo "Removing reads that failed the min overhang filter..."
awk -v fless=$lessfilename -v fmore=${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe 'BEGIN{i=0; k=0}{
    if ($0~/\-1\t\-1/){
         print $0 > fless;
         i++
    } 
    else{
         print $0 > fmore;
         k++;
    }
    }
END{
    printf ("Filtered out %f reads with <60% overlap with a single digestion fragment\n", i/(i+k));
}' ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe

if [ "$nodelete" != "nodelete" ]; then
        rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag.bedpe
fi

echo "Adding frag length and signed distance from bait; removing self-ligation fragments (if any; not expected with HiCUP input)..."
perl -ne '{ 
    chomp $_; 
    my @a = split /\t/, $_; 
    my $l = $a[16]-$a[15]; 
    my $d = "NA"; 
    if ($a[10] eq $a[14]){ 
           my $midB = int(($a[12]+$a[11])/2+0.5); 
           my $midR = int(($a[15]+$a[16])/2+0.5); 
           $d = $midR - $midB; 
    } 
    if ($d ne "0"){ 
           print "$_\t$l\t$d\n" 
    }
}' ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe > ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe

if [ "$nodelete" != "nodelete" ]; then
        rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06.bedpe
fi


echo "Pooling read pairs..."
echo "baitID	otherEndID	N	otherEndLen	distSign" > ${samplename}/${bamname}.chinput
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
}' ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe | sort -k1,1 -k2,2n -T ${samplename} >> ${samplename}/${bamname}.chinput

if [ "$nodelete" != "nodelete" ]; then
	rm ${samplename}/${bamname}_mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe
fi

echo "Done! The file to be used for Chicago R package input is ${samplename}/${bamname}.chinput"

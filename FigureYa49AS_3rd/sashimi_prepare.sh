#!/bin/bash

set -e
set -u
set -o pipefail

if [ $# -lt 4 ]
then
   echo -e "Usage: sashimi_prepare.sh -p sashimi -t 10 ref.fa isoform.fa read_1.fq read_2.fq"
   echo -e "\t-t: threads"
   echo -e "\t-p: output prefix"
   exit 1
fi

# check the software
for soft in gmap hisat2 samtools gffread
do
    if [ ! $(which $soft) ]
    then
        echo "please install $soft firstly" >&2
        exit 1
    fi
done

# set the parameter
THREADS=10
PREFIX=sashimi
while getopts t:p: opt
do
    case "$opt" in
    t) THREADS="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    *) echo "Unkown option: $opt";;
    esac
done

shift $[ $OPTIND -1 ]

REF=$1
ISOFORM=$2
FQ1=$3
FQ2=$4

mkdir -p index

# get the gff3 from isoform
## Build Index
if [ ! -f index/${PREFIX}/${PREFIX}.version ]; then
    gmap_build -D index -d $PREFIX $REF
fi

##Align
if [ ! -f ${PREFIX}.gff3 ]; then
    gmap -t $THREADS -D index  -d $PREFIX -f gff3_gene $ISOFORM > ${PREFIX}.gff3
fi

# covert gff3 to gtf for visualization
if [ ! -f ${PREFIX}.gtf ]; then
    gffread ${PREFIX}.gff3 -T -o ${PREFIX}_tmp.gtf
    sed -n '/exon/ p' ${PREFIX}_tmp.gtf > ${PREFIX}.gtf
    rm -f ${PREFIX}_tmp.gtf
fi

# Get the BAM
## Build Index
if [ ! -f index/${PREFIX}.1.ht2 ]; then
    hisat2-build -p $THREADS $REF index/${PREFIX}
fi

## Align
if [ ! -f align.done ]; then
    hisat2 -p $THREADS -x index/${PREFIX} -1 $FQ1 -2 $FQ2  | samtools sort -@ $THREADS > ${PREFIX}_sort.bam && \
    touch align.done
fi

## Filter
if [ ! -f filter.done ]; then
    samtools view -b -@ $THREADS -q 30 ${PREFIX}_sort.bam > ${PREFIX}_flt.bam && touch filter.done
fi

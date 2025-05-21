#! /bin/bash

#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH -t 7-0:00
#SBATCH -o spaceranger.out
#SBATCH -e spaceranger.err

while [[ $# -gt 0 ]];do
        case $1 in
                --sample) sample=$2; shift;;
                --slide) slide=$2; shift;;
                --area) area=$2; shift;;
                --thread) thread=$2; shift;;
                -*) shellcmd=$2;;
        esac
        shift
done

#if [[ ! -n $sample ]]; then sample="CRC_16"; fi
#if [[ ! -n $slide ]]; then slide="V10S15-020"; fi
#if [[ ! -n $area ]]; then area="D1"; fi
if [[ ! -n $thread ]]; then thread=8; fi

module load spaceranger

cmd="spaceranger count"
cmd="$cmd --id $sample --fastqs InputData/$sample --image $(ls InputData/$sample/*.tiff)"
cmd="$cmd --slide $slide --area $area --slidefile InputData/$sample/$slide.gpr"
cmd="$cmd --transcriptome Reference/refdata-gex-GRCh38-2020-A"
cmd="$cmd --localcores=$thread --localmem=$(($thread*10))"
echo $cmd; eval $cmd

mv $sample Results/
rm __$sample.mro

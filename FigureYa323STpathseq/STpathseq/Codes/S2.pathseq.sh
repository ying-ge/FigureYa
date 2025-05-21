#! /bin/bash

#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH -t 7-0:00
#SBATCH -o pathseq.out
#SBATCH -e pathseq.err

while [[ $# -gt 0 ]];do
        case $1 in
                --sample) sample=$2; shift;;
                -*) shellcmd=$2;;
        esac
        shift
done

# if [[ ! -n $sample ]]; then sample="CRC_16"; fi

gatk="/shared/home/lux/.conda/envs/gatk/bin/java -Xmx100G -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
cmd="$gatk PathSeqPipelineSpark"
cmd="$cmd --input Results/$sample/outs/possorted_genome_bam.bam"
cmd="$cmd --output Results/$sample/outs/$sample.pathseq.complete.bam"
cmd="$cmd --scores-output Results/$sample/outs/$sample.pathseq.complete.csv"
cmd="$cmd --filter-bwa-image Reference/pathseq_host/pathseq_host.fa.img"
cmd="$cmd --microbe-dict Reference/pathseq_microbe/pathseq_microbe.dict"
cmd="$cmd --microbe-bwa-image Reference/pathseq_microbe/pathseq_microbe.fa.img"
cmd="$cmd --taxonomy-file Reference/pathseq_microbe/pathseq_taxonomy.db"
cmd="$cmd --kmer-file Reference/pathseq_host/pathseq_host.bfi"
cmd="$cmd --min-clipped-read-length 60 --min-score-identity .7"
cmd="$cmd --is-host-aligned false --filter-duplicates false"

echo $cmd; eval $cmd

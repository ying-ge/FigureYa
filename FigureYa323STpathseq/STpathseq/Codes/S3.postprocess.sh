#! /bin/bash

#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH -t 7-0:00
#SBATCH -o post.out
#SBATCH -e post.err

module load r/4.1.1

while [[ $# -gt 0 ]];do
        case $1 in
                --sample) sample=$2; shift;;
                -*) shellcmd=$2;;
        esac
        shift
done

# if [[ ! -n $sample ]]; then sample="CRC_16"; fi

cmd="python Codes/UMI_annotator.py --sample $sample"
echo $cmd; eval $cmd

cmd="python Codes/validate_and_count.py --sample $sample"
# echo $cmd; eval $cmd

cmd="Rscript Codes/SeuratAnalysis.R --samplename $sample"
echo $cmd; eval $cmd

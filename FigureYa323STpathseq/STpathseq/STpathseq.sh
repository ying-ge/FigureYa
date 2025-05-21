#! /bin/bash

#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH -t 7-0:00
#SBATCH -o STpathseq.out
#SBATCH -e STpathseq.err

module load spaceranger
module load r/4.1.1

for i in $(cat InputData/sample_sheet.csv|grep -v \#|wc -l|xargs seq)
do
{
  tmp=$(cat InputData/sample_sheet.csv|grep -v \#|sed -n ${i}p)
  sample=$(echo $tmp|cut -d, -f1)    
  slide=$(echo $tmp|cut -d, -f2)
  area=$(echo $tmp|cut -d, -f3)

  echo
  echo "========= sample: $sample ==========="
  echo "=========  slide: $slide  ==========="
  echo "=========   area: $area   ==========="

    bash Codes/S1.spaceranger.sh --sample $sample --slide $slide --area $area
    bash Codes/S2.pathseq.sh --sample $sample
    bash Codes/S3.postprocess.sh --sample $sample

  echo "======================================"
  echo
  echo
}
done

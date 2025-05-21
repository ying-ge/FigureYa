cd ../data/

cellranger count --id=old1 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/O1/ \
--sample=O1_Sample_SaS_CB_002 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64

cellranger count --id=old2 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/O2/ \
--sample=O2_SaS_CB_005 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64

cellranger count --id=old3 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/O3/ \
--sample=O3_SaS_CB_011 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64

cellranger count --id=young1 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/Y1/ \
--sample=Y1_Sample_SaS_CB_001 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64

cellranger count --id=young2 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/Y2/ \
--sample=Y2_SaS_CB_003 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64

cellranger count --id=young3 \
--transcriptome=../ref/gencode_vM12 \
--fastqs=fastq/Y3/ \
--sample=Y3_SaS_CB_004 \
--localcores=16 \
--expect-cells=10000 \
--localmem=64
#!bin/bash
ARACNE_BIN=/ifs/scratch/c2b2/ac_lab/lv2395/ARACNe/aracne-ap.jar
ARACNE_DIR='/ifs/scratch/c2b2/ac_lab/lv2395/ARACNe'
BOOTSTRAP_NUM=200
P_THRESH=1E-8

## initial setting of cleanup, reg file list
CLEANUP=true
REG_FILES=()
## process arguments
while (( "$#" )); do
	case "$1" in
		-n|--run-id)
			RUN_ID=$2
			shift 2
			;;
		-b|--base-dir)
			BASE_DIR=$2
			shift 2
			;;
		-i|--input_file)
			IN_FILE=$2
			shift 2
			;;
		-r|--reg-set)
			REG_FILES=(${REG_FILES[@]} $2)
			shift 2
			;;
		-k)
			CLEANUP=false
			shift 1
			;;
		-*|--*-)
			echo 'Error: Unsupported flag $1'
			exit 1
			;;
	esac
done
## create regulator sub-directories
REG_SETS=()
for (( c = 0; c < ${#REG_FILES[@]}; c++ ))
do
	# parse string
	FILE=${REG_FILES[$c]}
	REG=${FILE##*/}
	REG=${REG%%.*}
	# add to list and make directory
	mkdir -p ${BASE_DIR}/${RUN_ID}_${REG}
	REG_SETS=(${REG_SETS[@]} ${REG})
done
## make log folder
LOG_DIR=${BASE_DIR}/${RUN_ID}_log-files
mkdir -p ${LOG_DIR}

## Prep data from .rds to ARACNe Table
A_TABLE=${BASE_DIR}${RUN_ID}_ARACNe-table.tsv
qsub -N adp_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=:30: -l h_rt=:30: \
-o ${LOG_DIR}/${RUN_ID}_adp-o.txt -e ${LOG_DIR}/${RUN_ID}_adp-e.txt -b y \
/nfs/apps/R/3.5.1/bin/Rscript ${ARACNE_DIR}/aracne_data-prep.r --rds=${IN_FILE} --out=${A_TABLE}

## Calculate MI, generate boostraps, and consolidate for each reg set
HOLD_LIST=''
for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
do
	## Set iteration variables
	REG_FILE=${REG_FILES[c]}
	REG=${REG_SETS[c]}
	WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
	RUN_NAME=${RUN_ID}_${REG}
	HOLD_LIST="$HOLD_LIST,aCon_${RUN_NAME}"
	## Threshold calculation
	qsub -N aThresh_${RUN_NAME} -hold_jid adp_${RUN_ID} -wd ${WORK_DIR} -l mem=16G -l s_rt=2:0:0 -l h_rt=2:0:0 -b y \
	java -Xmx5G -jar ${ARACNE_BIN} -e ${A_TABLE} -o ${WORK_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed 666 --calculateThreshold
	## Bootstrap generation
	for (( i = 1; i <= $BOOTSTRAP_NUM; i++ ))
	do
		qsub -N abs_${RUN_NAME} -hold_jid aThresh_${RUN_NAME} -wd ${WORK_DIR} -l mem=16G -l s_rt=4:0:0 -l h_rt=4:0:0 -b y \
		java -Xmx5G -jar ${ARACNE_BIN} -e ${A_TABLE} -o ${WORK_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed $i
	done
	## Consolidation
	qsub -N aCon_${RUN_NAME} -hold_jid abs_${RUN_NAME} -wd ${WORK_DIR} -l mem=16G -l s_rt=8:0:0 -l h_rt=8:0:0 -b y \
	/nfs/apps/R/3.5.1/bin/Rscript ${ARACNE_DIR}/aracne_consolidate.r ${WORK_DIR} ${A_TABLE} ${REG_FILE} bonferroni 0.01
done 

## Create Hold List
CAT_LIST=''
for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
do
	REG=${REG_SETS[c]}
	WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
	CAT_LIST="$CAT_LIST ${WORK_DIR}/finalNetwork_4col.tsv"
done
HOLD_LIST=${HOLD_LIST#?}
## Merge
qsub -N aMerge_${RUN_ID} -hold_jid ${HOLD_LIST} -wd ${BASE_DIR} -l mem=16G -l s_rt=2:0:0 -l h_rt=2:0:0 \
-o ${LOG_DIR}/${RUN_ID}_amerge-o.txt -e ${LOG_DIR}/${RUN_ID}_amerge-e.txt -b y \
bash ${ARACNE_DIR}/aracne_merge.sh ${CAT_LIST} ${BASE_DIR}${RUN_ID}_finalNet-merged.tsv
## Reg Process
qsub -N arp_${RUN_ID} -hold_jid aMerge_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 \
-o ${LOG_DIR}/${RUN_ID}_arp-o.txt -e ${LOG_DIR}/${RUN_ID}_arp-e.txt -b y \
/nfs/apps/R/3.5.1/bin/Rscript ${ARACNE_DIR}/aracne_regProc.r --a_file=${BASE_DIR}${RUN_ID}_finalNet-merged.tsv --exp_file=${A_TABLE} --out_dir=${BASE_DIR} --out_name=${RUN_ID}
## Cleanup
if ${CLEANUP}
then
	qsub -N aClean_${RUN_ID} -hold_jid arp_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=1:0:0 -l h_rt=1:0:0 \
	-o ${LOG_DIR}/${RUN_ID}_acu-o.txt -e ${LOG_DIR}/${RUN_ID}_acu-e.txt -b y \
	bash ${ARACNE_DIR}/aracne_cleanup.sh ${BASE_DIR} ${RUN_ID} ${REG_SETS[@]}
fi

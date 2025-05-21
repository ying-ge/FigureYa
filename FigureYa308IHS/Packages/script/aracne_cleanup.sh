#!bin/bash

## get arguments
BASE_DIR=$1
RUN_ID=$2
REG_SETS=( ${@:3} )

# ## make file for meta data / folder for .tsvs
META_FILE=${BASE_DIR}/${RUN_ID}_metaData.csv
echo 'Regulator,Bootstraps,Interactions' >> ${META_FILE}
mkdir -p ${BASE_DIR}/${RUN_ID}_net-tsvs
## clean up each regulator
for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
do
	REG=${REG_SETS[c]}
	NET_FILE=${BASE_DIR}${RUN_ID}_${REG}/finalNetwork_4col.tsv
	## write metadata
	boot_num=$(ls ${BASE_DIR}${RUN_ID}_${REG}/bootstrapNetwork_* | wc -l)
	edge_num=( $(wc -l ${NET_FILE}) )
	echo ${REG},${boot_num},${edge_num} >> ${META_FILE}
	## move tsv
	mv ${NET_FILE} ${BASE_DIR}/${RUN_ID}_net-tsvs/${RUN_ID}_${REG}_4col.tsv
	## delete directory
	rm -r ${BASE_DIR}/${RUN_ID}_${REG}
done

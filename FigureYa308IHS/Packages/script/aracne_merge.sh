#!bin/bash
NUM_ARGS=$#
CAT=''
for (( c = 1 ; c < NUM_ARGS ; c++ ))
do
	CAT="$CAT ${!c}"
done
cat ${CAT} > ${!NUM_ARGS}
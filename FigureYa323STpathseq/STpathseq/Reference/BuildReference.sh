#! /bin/bash

#SBATCH -t 7-0:00
#SBATCH -o reference.out
#SBATCH -e reference.err

# if [ ! -d Reference ];then mkdir Reference; fi

md5sum -c md5sum.txt > log

if [ $(cat log|grep pathseq|grep -v OK|wc -l) -gt 0 ]; then
	prefix="wget -r ftp://gsapubftp-anonymous:@"
	for file in $(cat log|grep pathseq|grep -v OK|cut -d: -f1)
	do
		cmd="$prefix$file"; echo $cmd; eval $cmd
	done
fi

if [ $(cat log|grep refdata|grep -v OK|wc -l) -gt 0 ];then
	cmd="wget https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"
	echo $cmd; eval $cmd
fi

echo uncompress the reference

if [ ! -d pathseq_host ]; then mkdir pathseq_host; fi
if [ ! -d pathseq_microbe ];then mkdir pathseq_microbe; fi

tar -zxkvf ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz -C pathseq_host
tar -zxkvf ftp.broadinstitute.org/bundle/pathseq/pathseq_taxonomy.tar.gz -C pathseq_microbe
tar -zxkvf ftp.broadinstitute.org/bundle/pathseq/pathseq_microbe.tar.gz -C pathseq_microbe
tar -zxkvf refdata-gex-GRCh38-2020-A.tar.gz #-C Reference
#mv Reference/pathseq_host.bfi Reference/pathseq_host/

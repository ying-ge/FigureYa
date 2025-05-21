# target="../data/processed"
# mkdir -p $target
# cat processed.txt | sort | uniq | xargs -I {} wget -P $target {}

# cd $target
# find -name "*.zip" | xargs -I {} unzip {}
# cd -

# download fastq
base=`pwd`
target="../data/raw"
mkdir -p $target
cd $target

sed '1d' $base/fastq.txt | awk '{print "mkdir -p",$1}' | bash
sed '1d' $base/fastq.txt | awk '{print "wget -P",$1,$2,"\n","wget -P",$1,$3,"\n","wget -P",$1,$4}' | python $base/run_parallel.py 9
# sed '1d' $base/fastq.txt | awk '{print "mkdir",$1,"\n","wget -P",$1,$2,"\n","wget -P",$1,$3,"\n","wget -P",$1,$4}' | bash

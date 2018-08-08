# @Author: Francis Blokzijl
# @Date: 15 Mar 2017
# @Description: Determine surveyed area of autosomal mouse chromosomes
# bed files are already intersected between control sample and subject sample
# 1. cat bed file
# 2. select autosomal mouse chromosomes
# 3. add up ranges

# change dir
cd ~/surfdrive/Shared/ERCC1/Data/callable_genome/bed
# list all files in dir
FILES=`ls`

for f in $FILES
do
	echo $f >> ~/surfdrive/Shared/ERCC1/Data/callable_genome/determine_callable/callable_autosomal.txt 
	cat $f | awk '{if($1 >= 1 && $1 <= 19) {print $0}}' | awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' >> ~/surfdrive/Shared/ERCC1/Data/callable_genome/determine_callable/callable_autosomal.txt
done 
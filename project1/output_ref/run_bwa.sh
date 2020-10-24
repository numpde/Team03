# RA, 2020-10-24

I=../input/data
O=data
mkdir -p ${O}
bwa index ${I}/genome.chr22.fa.gz -p ${O}/large
bwa mem ${O}/large ${I}/*5xCov1.fq.gz ${I}/*5xCov2.fq.gz > ${O}/large.sam

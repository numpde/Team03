# RA, 2020-09-28
I=../input/data_small
O=data_small
mkdir -p ${O}
bwa index ${I}/genome.chr22.5K.fa -p ${O}/tiny
bwa mem ${O}/tiny ${I}/*30xCov1.fq ${I}/*30xCov2.fq > ${O}/tiny.sam


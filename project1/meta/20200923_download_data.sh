#!/bin/env bash

# Download the /Comp Biomed, Project 1/ input files
# RA, 2020-09-23
# https://gist.github.com/numpde/dcd97395e9d12b9fb2e5113f688f5339

ORIG=$(pwd)

wget -N https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project1/project1_description.pdf

cd "${ORIG}"
DIR=data_small
mkdir -p ${DIR}
cd ${DIR}
BASE=https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project1/${DIR}
wget -N ${BASE}/genome.chr22.5K.fa
wget -N ${BASE}/output_tiny_30xCov1.fq
wget -N ${BASE}/output_tiny_30xCov2.fq
wget -N ${BASE}/output_tiny_30xCov.sam


cd "${ORIG}"
DIR=data
mkdir -p ${DIR}
cd ${DIR}
BASE=https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project1/${DIR}
wget -N ${BASE}/genome.chr22.fa.gz
wget -N ${BASE}/output_5xCov1.fq.gz
wget -N ${BASE}/output_5xCov2.fq.gz
wget -N ${BASE}/output_10xCov1.fq.gz
wget -N ${BASE}/output_10xCov2.fq.gz
wget -N ${BASE}/output_30xCov1.fq.gz
wget -N ${BASE}/output_30xCov2.fq.gz

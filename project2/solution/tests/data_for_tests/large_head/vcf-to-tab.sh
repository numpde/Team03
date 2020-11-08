#!/usr/bin/env bash

# RA, 2020-11-08

# sudo apt install vcftools

vcf-to-tab  < case_processed.vcf  > case_processed.vcf.tab
vcf-to-tab  < control.vcf         > control.vcf.tab

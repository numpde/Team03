# WARNING! Don't run this script, it will download a 15GB and expand to 110GB

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz

gunzip GRCh37_latest_dbSNP_all.vcf.gz

grep '^NC_000017.10' GRCh37_latest_dbSNP_all.vcf > temp.vcf

head -n38 GRCh37_latest_dbSNP_all.vcf > GRCh37_latest_dbSNP_chrom17.vcf
cat temp.vcf >> GRCh37_latest_dbSNP_chrom17.vcf
rm temp.vcf
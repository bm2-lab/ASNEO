#!/usr/bin/env sh
if [ "$#" -ne 1 ];then
	echo -e "Test ASNEO\nUsage: run_ASNEO.sh reference_genome(hg19 or GRCh37)\n"
	exit 1
fi

ref_genome=$1   #~/database/hg19/hg19.fa
python ../ASNEO.py -j SRR2660032.SJ.out.tab -a HLA-C05:01,HLA-B18:01,HLA-C07:02,HLA-A03:01,HLA-A30:02,HLA-B07:02 -g $ref_genome 

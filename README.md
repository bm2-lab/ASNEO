# _ASNEO: Alternative Splicing NEOantigens_

## Introduction
ASNEO is a computational pipeline to identify personalized tumor neoantigens derived from alternative splicing.

## Requirement
### software
* python
* bedtools
* pepmatch_db_x86_64 (a script in MuPeXI)
* netMHCpan
* netCTLpan

### python package
* pandas
* sj2psi
* pysam
* biopython

## Usage
### input
1. Reference genome file: `hg19.fa` or `GRCh37.fa`
2. `STAR` mapped file: `SJ.out.tab` and `indexed bam`(optional)

### example
1. git clone https://github.com/zzb23/ASNEO.git
2. cd ASNEO/src && tar -xzvf software.tar.gz 
3. cd ../example && bash run_ASNEO.sh hg19.fa 
> `hg19.fa` should change to your own reference genome file (must be hg19 or GRCh37)

For detailed usage information, please refer to the [ASNEO User Manual](/doc/ASNEO_User_Manual.md)

## Citation
TODO

## Contact
qiliu@tongji.edu.cn

Tongji University, Shanghai, China

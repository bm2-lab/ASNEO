# _ASNEO: identification of Alternative Splicing NEOantigens_

## Introduction
ASNEO is a computational pipeline to identify personalized tumor neoantigens derived from alternative splicing.

## Requirement
### software
* python
* bedtools
* pepmatch_db_x86_64 (a script in MuPeXI)
* netMHCpan
* netCTLpan
> The software are packed in `src/software.tar.gz`

### python package
* pandas
* sj2psi
* pysam
* biopython
* xgboost

## Usage

### Parameters

* ##### Required

  * -j, --junc

	Junction file, the default format is STAR result SJ.out.tab.  You can also specify your own format by columns names in your junction file and setting --columns. 

  * -g, --genome

	Reference genome file, must be hg19 or GRCh37 .

* ##### Optional

  * -l, --length

	Peptide length [8-11], default "9". if multiple, seperated with ",", such as 8,9,10,11.

  * -a, --allele

	HLA allele, default "HLA-A02:01". Multiple separated with ",". You can use OptiType or Seq2HLA to detect it with RNA-Seq.

  * -p, --process

	Multiple processes, default "10".

  * -o, --outdir

	Output directory, default "result".

  * --reads

	Cutoff of junction reads number, default "10". Only junctions whose unique reads number greater than it will remain.

  * --psi

	Cutoff of psi value(0-1), default 0.1.  [See PSI means here](https://github.com/olgabot/sj2psi).

  * --bind

	Cutoff of bind rank,  default 2. 0.5 was high binding peptides rank threshold(SB), 2 was low binding peptides rank threshold(WB).

  * --bam

	Bam file, default None. Our software applies a simple calculation of RPKM to filter neoantigens. If you want to use it, specify the indexed bam file. You can also filter the neoantigens by yourself after process done with your own express calculation. 

  * --rpkm

	Cutoff of RPKM value, default 1.0. Only useful when you set --bam.

  * --columns

	If junction file contains columns names. If set, at least contain "chrom, intron_start, intron_stop, unique_junction_reads". Note: the intron coordinates must be 1-based and columns must be separated by Tab.


### Example

1. git clone https://github.com/bm2-lab/ASNEO.git

2. cd ASNEO/src && tar -xzvf software.tar.gz

3. cd ../test && bash run_ASNEO.sh hg19.fa
> `hg19.fa` should change to your own reference genome file (must be hg19 or GRCh37)

### Output

The output file "putative_neo.txt" contains all putative neoantigens information.

| Column | Description |
| - | -: |
| HLA | HLA type |
| mismatch | The base number different between mtpep and wtpep |
| MTpep | Alternative splicing derived neopeptide |
| MTpep_score | netMHCpan predicted score of mtpep |
| MTpep_aff | netMHCpan predicted affinity of mtpep |
| MTpep_rank | netMHCpan predicted bind rank(%) of mtpep |
| MTpep_comb | netCTLpan predicted comb score of mtpep, combined MHC score, cleavage score and TAP score |
| WTpep | pepmatch_db_x86_64 extracted normal peptide |
| WTpep_score |  netMHCpan predicted score of wtpep |
| WTpep_aff | netMHCpan predicted affinity of wtpep |
| WTpep_rank | netMHCpan predicted bind rank(%) of wtpep |
| WTpep_comb | netCTLpan predicted comb score of wtpep, combined MHC score, cleavage score and TAP score |
| Hydro_Model | Neopeptide immunity measurement based on amino acid hydrophobicity |
| R | T cell recognition score calculated based on TCR cross-reactivity |
| Score | Immunogenicity score for neoantigens |
| Junc_info | The source junctions' information of the neoantigen, the format like this:<br>transcript_id:protein_id:gene_symbol@JUNC:chrom_start_stop@RPKM:rpkm_value(strand). <br>multiple junctions are separated by "\|". <br>Note: start,stop are 0-based intron start-stop coordinate; rpkm_value=-1 means we didn't calculate it  |


## Citation
ASNEO: identification of personalized Alternative Splicing based NEOantigens with RNA-seq, submitted, 2019

## Contact
zhangzb@tongji.edu.cn or qiliu@tongji.edu.cn

Tongji University, Shanghai, China

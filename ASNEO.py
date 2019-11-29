#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function
import os, logging, sj2psi, sys, subprocess, shutil, pysam
import pickle, tempfile, multiprocessing, warnings
from subprocess import PIPE
import pandas as pd
from functools import partial
from math import log, exp
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from argparse import ArgumentParser
import sklearn
from xgboost import XGBClassifier

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]:  %(message)s')


def InputParser():
    parser = ArgumentParser()
    parser.add_argument('-j', '--junc', required=True, help='Input junction file')
    parser.add_argument('-g', '--genome', required=True, default='hg19.fa', help='Referece genome file')
    parser.add_argument('-l', '--length', default='9', help='Peptide length [8-11] (multiple length with ,)')
    parser.add_argument('-a', '--allele', default='HLA-A02:01', help='HLA allele', dest='allele')
    parser.add_argument('-p', '--process', default=10, help='Multiprocess number to use')
    parser.add_argument('-o', '--outdir', default='result', help='Output directory')
    parser.add_argument('--reads', default=10, help='Cutoff of junction reads number, default 10')
    parser.add_argument('--psi', default=0.1, help='Cutoff of psi value(0-1), default 0.1')
    parser.add_argument('--bind', default=2, help='Cutoff of bind rank,SB:2,WB:0.5,default 2')
    parser.add_argument('--bam', default=None, help='Bam file if you want to calculate rpkm')
    parser.add_argument('--rpkm', default=1, help='Cutoff of transcript rpkm value, default 1')
    parser.add_argument('--columns', action='store_true', help='If junction file contain columns name,'
                             'If set, must contain "chrom, intron_start, intron_stop, unique_junction_reads"')
    args = parser.parse_args()
    return args


def FormatJunc(junc_file, columns):
    if not columns:
        COLUMN_NAMES = ('chrom', 'intron_start', 'intron_stop', 'strand', 'intron_motif',
                        'annotated', 'unique_junction_reads', 'multimap_junction_reads', 'max_overhang')
        try:
            sj = pd.read_table(junc_file, header=None, names=COLUMN_NAMES, sep='\s+')
        except:
            logging.error('Load junction file failed, please check the format!')
            sys.exit(1)
    else:
        try:
            sj = pd.read_table(junc_file, header=0, sep='\s+')
            sj['multimap_junction_reads'] = 0
        except:
            logging.error('Load junction file failed, please check the format!')
            sys.exit(1)

    sj['first_bp_intron'] = sj.intron_start
    sj['last_bp_intron'] = sj.intron_stop
    all_chroms = ['chr' + str(i + 1) for i in range(22)] + ['chrX', 'chrY']
    sj = sj[sj.chrom.astype(str).isin(all_chroms)]
    return sj


def FilterSJ(sj, reads, psi):
    # filter psi
    sj = sj2psi.get_psis(sj, min_unique=int(reads), min_multimap=0)
    sj = sj[(sj.psi5 > float(psi)) & (sj.psi3 > float(psi))]
    # filter count
    sj = sj[sj['unique_junction_reads'] > int(reads)]
    # filter gtex normal junctions
    gtex = pd.read_table(path['gtex'], header=None)
    sj.index = sj['chrom'].str.strip('chr') + '_' + sj['intron_start'].astype(str) + '_' \
               + sj['intron_stop'].astype(str)
    sj = sj[~sj.index.isin(gtex.iloc[:, 0])]
    return sj


def CalRPKM(iso, bam):
    bam = pysam.AlignmentFile(bam, 'rb')
    all_chroms = ['chr' + str(i + 1) for i in range(22)] + ['chrX', 'chrY']
    bam_size = sum([stat.mapped for stat in bam.get_index_statistics() if stat.contig in all_chroms])
    count = bam.count(contig=iso.chrom, start=iso.txStart + 1, stop=iso.txEnd)
    rpkm = count / sum(map(int, iso.exonLens.strip(',').split(','))) / bam_size * 1E9
    return round(rpkm, 2)


def Jun2Iso(junc, proteome, bam, rpkm_value):
    chrom, intron_start, intron_end = junc['chrom'], junc['intron_start'] - 1, junc['intron_stop']  # zero-based
    isoforms = proteome[(proteome.chrom == chrom) & (proteome.txStart < intron_start) & (proteome.txEnd > intron_end)]
    # use tag to define the junction, 0: annotated junction; 1: novel junction;
    # 2: novel donor site or novel acceptor site; 3: both novel site
    tag_iso = {}
    novel_isoforms = set()
    if not isoforms.empty:  # junction in transcript region
        for _, iso in isoforms.iterrows():  # for each iso, search the junction
            exon_starts = list(map(int, iso.exonStarts.strip(',').split(',')))
            exon_lens = list(map(int, iso.exonLens.strip(',').split(',')))
            exon_ends = [exon_start + exon_len for exon_start, exon_len in zip(exon_starts, exon_lens)]
            novel_exon_end = intron_start - iso.txStart
            novel_exon_start = intron_end - iso.txStart
            # define the iso tag
            u_index = exon_ends.index(novel_exon_end) if novel_exon_end in exon_ends else -1  # upstream exon index
            d_index = exon_starts.index(
                novel_exon_start) if novel_exon_start in exon_starts else -1  # downstream exon index
            if u_index > -1 and d_index > -1:  # two site exists
                tag = 0 if d_index - u_index == 1 else 1
            elif u_index == -1 and d_index > -1:
                tag = 2
            elif u_index > -1 and d_index == -1:
                tag = 2
            else:
                tag = 3
            tag_iso.setdefault(tag, []).append(iso)
    if tag_iso and 0 not in tag_iso:
        min_intron = 50  # minimum intron length for novel junctions
        min_exon = 2  # minmum exon length changed by novel exon in both upstream and downstream
        max_exon = 250  # maximum exon length expanded by novel exon in both upstream and downstream
        tag = sorted(tag_iso.keys())[0]  # extract lowest tag
        for iso in tag_iso[tag]:
            if intron_start > iso.cdsStart and intron_end < iso.cdsEnd:  # junction not in protein coding region, abort!!
                # shift transcript exon coordinate to cds exon coordinate
                abs_exon_starts = [exon_start + iso.txStart for exon_start in
                                   list(map(int, iso.exonStarts.strip(',').split(',')))]
                exon_lens = list(map(int, iso.exonLens.strip(',').split(',')))
                abs_exon_ends = [exon_start + exon_len for exon_start, exon_len in zip(abs_exon_starts, exon_lens)]

                u_index = [i for i, coord in enumerate(abs_exon_starts) if int(coord) < intron_start][-1]
                d_index = [i for i, coord in enumerate(abs_exon_ends) if int(coord) > intron_end][0]
                judge = [intron_end - intron_start, intron_start - abs_exon_ends[u_index],
                         abs_exon_starts[d_index] - intron_end,
                         intron_start - abs_exon_starts[u_index], abs_exon_ends[d_index] - intron_end]
                if judge[0] > min_intron and judge[1] < max_exon and judge[2] < max_exon \
                        and abs(judge[1]) > min_exon and abs(judge[2]) > min_exon \
                        and judge[3] > min_exon and judge[4] > min_exon:
                    novel_exon_starts = abs_exon_starts[0:u_index + 1] + [intron_end] + abs_exon_starts[d_index + 1:]
                    novel_exon_ends = abs_exon_ends[0:u_index] + [intron_start] + abs_exon_ends[d_index:]
                    cds_exon_starts = [0] + [start - iso.cdsStart for start in novel_exon_starts if
                                             iso.cdsStart < start < iso.cdsEnd]
                    cds_exon_ends = [end - iso.cdsStart for end in novel_exon_ends if
                                     iso.cdsStart < end < iso.cdsEnd] + [iso.cdsEnd - iso.cdsStart]
                    cds_exon_lens = [exon_end - exon_start for exon_end, exon_start in
                                     zip(cds_exon_ends, cds_exon_starts)]
                    rpkm = CalRPKM(iso,bam) if bam else -1
                    if rpkm == -1 or rpkm > rpkm_value:
                        name = '%s:%s:%s@JUNC:%s_%s_%s@RPKM:%s' % \
                               (iso.isoform, iso.protein, iso.gene, chrom, intron_start, intron_end, rpkm)
                        novel_iso = [iso.chrom, iso.cdsStart, iso.cdsEnd, name, iso.protein, iso.strand, iso.txStart,
                                     iso.txEnd, iso.gene,
                                     len(cds_exon_starts), ','.join(map(str, cds_exon_lens)),
                                     ','.join(map(str, cds_exon_starts))]
                        novel_iso = '\t'.join(map(str, novel_iso))
                        novel_isoforms.add(novel_iso)
    return novel_isoforms


def MultiIters(its, func, process):
    num = 0
    for i in range(len(its)//200 + 1):
        pool = multiprocessing.Pool(processes=int(process))
        bits = (it[1] for it in its.iloc[i*200:(i+1)*200,].iterrows())
        for out in pool.imap(func, bits):
            yield out
            num += 1
            sys.stdout.write('Done %d/%d\r' % (num, len(its)))
        pool.close()
        pool.join()


def ProcessSJ(junc_file, columns, reads, psi, bam, rpkm, process):
    logging.info("Begin Junction to Isoform Process ...")
    sj = FormatJunc(junc_file=junc_file, columns=columns)
    logging.info("Origin junctions number is: %s", len(sj))
    sj = FilterSJ(sj=sj,reads=reads,psi=psi)
    logging.info("After filter, Remain junctions number is: %s", len(sj))

    logging.info("Generate novel isoforms")
    proteome = pd.read_table(path['refseq_bed'], header=None, comment='#', names=('chrom', 'txStart',
    'txEnd', 'isoform', 'protein', 'strand','cdsStart', 'cdsEnd', 'gene', 'exonNum','exonLens', 'exonStarts'))
    isoforms, count = set(), 0
    novel_isoforms = MultiIters(its=sj, func=partial(Jun2Iso, proteome=proteome, bam=bam, rpkm_value=rpkm), process=process)
    for iso in novel_isoforms:
        if iso:
            isoforms = isoforms | iso
            count += 1
    logging.info("Novel isoforms number is %s" , len(isoforms))

    with open(path['iso_bed'], 'w') as f:
        f.write('\n'.join(isoforms) + '\n')


def Iso2Prot(genome, isobed, isofa, protfa):
    logging.debug("Call bedtools to extract DNA sequence")
    cmd = "%s getfasta -s -name -split -fi %s -bed %s -fo %s" % (path['bedtools'], genome, isobed, isofa)
    subprocess.call(cmd, shell=True)
    seqs = SeqIO.parse(isofa, 'fasta')
    with open(protfa, 'w') as f:
        for sequence in seqs:
            pep = sequence.seq[0:len(sequence) // 3 * 3].translate()
            pep = pep.split('*')[0]
            if len(pep.strip()) > 30:
                f.write('>%s\n%s\n' % (sequence.id, pep))


def GenKmerPep(protfa, peplen, save=None):
    sequences = SeqIO.parse(protfa, 'fasta')
    peps = set()
    for sequence in sequences:
        sequence = str(sequence.seq).strip()
        for i in range(len(sequence) - peplen + 1):
            peps.add(sequence[i:i + peplen])
    if save:
        with open(save, 'w') as f:
            f.write('\n'.join(peps))
    return peps


def RunMHCpan(allele, peptxt, affit):
    cmd = '%s -a %s -p %s -BA > %s' % (path['netMHCpan'], allele, peptxt, affit)
    subprocess.call(cmd, shell=True)
    with open(affit) as f:
        line = f.readlines()[-1].strip()
        if "Error" in line:
            logging.error("Some error in call netMHCpan: %s"%line)
            sys.exit(1)


def ParseAffit(protfa, affit, epit, bind):
    sequences = SeqIO.to_dict(SeqIO.parse(protfa, 'fasta'))
    with open(affit) as fin, open(epit, 'w') as fout:
        for line in fin:
            lines = line.strip().split()
            if lines and lines[0] == '1' and float(lines[13]) < float(bind):
                hla, mtpep = lines[1].replace('*', ''), lines[2]
                keys = [list(sequences.keys())[i] for i, seq in enumerate(list(sequences.values())) if mtpep in seq]
                fout.write('\t'.join([hla, mtpep] + lines[11:14] + ['|'.join(keys)]) + '\n')


def ProcessIsoform(genome, length, allele, bind):
    logging.info("Begin Isoform to Epitope Process ...")
    logging.info("Translate isoforms to protein")
    Iso2Prot(genome=genome, isobed=path['iso_bed'], isofa=path['iso_fa'], protfa=path['prot_fa'])
    logging.info("Cut protein to short peptide")
    remain_pep = set()
    peplens = [int(l.strip()) for l in str(length).strip().split(',')]
    for peplen in peplens:
        peps = GenKmerPep(protfa=path['prot_fa'], peplen=peplen)
        norm_peps = GenKmerPep(protfa=path['refseq_fa'], peplen=peplen, save=path['refpep_%s'%peplen])
        peps = peps - norm_peps
        logging.debug("After filter, remain %s mer peptide number is: %s", peplen, len(peps))
        remain_pep = remain_pep | peps
    logging.info("All short pep number is: %s", len(remain_pep))
    with open(path['pep_txt'], 'w') as f:
        f.write('\n'.join(remain_pep) + '\n')

    logging.info("Call netMHCpan")
    RunMHCpan(allele=allele, peptxt=path['pep_txt'], affit=path['affit'])
    ParseAffit(protfa=path['prot_fa'], affit=path['affit'], epit=path['epit'], bind=bind)


def gentmp(line,tmpdir):
    tmp = tempfile.mkstemp(dir=tmpdir, text=True)[1]
    with open(tmp, 'w') as fout:
        fout.write(line)
    return tmp


def GetWildpep(mut_pep):
    line = '{}\n'.format(mut_pep)
    tmp = gentmp(line,path['tdir'])
    cmd = '%s -thr 10 %s %s' % (path['pepmatch'], tmp, path['refpep_%s'%len(mut_pep)])
    p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE).stdout
    lines = [line.decode() for line in p if not line.startswith(b'#')][0].strip().split()
    os.remove(tmp)
    return lines[3], lines[5]


def GetCTLpan(allele, pep):
    line = '>pep\n{}\n'.format(pep)
    tmp = gentmp(line,path['tdir'])
    cmd = '%s -a %s -l %s -f %s' % (path['netCTLpan'], allele, len(pep), tmp)
    p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE).stdout
    lines = [line for line in p if pep in line.decode()][0].decode().strip().split()
    os.remove(tmp)
    return lines[3], lines[4], lines[5]


def GetMHCpan(allele, pep):
    line = '>pep\n{}\n'.format(pep)
    tmp = gentmp(line,path['tdir'])
    cmd = '%s -a %s -l %s -BA -f %s' % (path['netMHCpan'], allele, len(pep), tmp)
    p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE).stdout
    lines = [line for line in p if pep in line.decode()][0].decode().strip().split()
    os.remove(tmp)
    return lines[11], lines[12], lines[13]


def aligner(seq1, seq2):
    matrix = matlist.blosum62
    gap_open = -11
    gap_extend = -1
    aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
    return aln


def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x - ma), v))) + ma


def GetR(pep):
    iedb_seq = SeqIO.parse(path['iedb'], 'fasta')
    align_score, a, k = [], 26, 4.86936
    for iseq in iedb_seq:
        aln_score = aligner(pep, iseq.seq)
        if aln_score:
            localds_core = max([line[2] for line in aln_score])
            align_score.append(localds_core)
    bindingEnergies = list(map(lambda x: -k * (a - x), align_score))
    lZk = logSum(bindingEnergies + [0])
    lGb = logSum(bindingEnergies)
    R = exp(lGb - lZk)
    return R


def do_hydro_vector(pep):
    hydro_score = dict(R=-4.5, K=-3.9, N=-3.5, D=-3.5, Q=-3.5, E=-3.5, H=-3.2, P=-1.6, Y=-1.3, W=-0.9, S=-.8,
                       T=-0.7, G=-0.4, A=1.8, M=1.9, C=2.5, F=2.8, L=3.8, V=4.2, I=4.5)
    hydro_vector = []
    for i in pep:
        hydro_vector.append(hydro_score[i.upper()])
    return hydro_vector


def GetH(pep):
    hydro_vector = do_hydro_vector(pep)
    if len(hydro_vector) == 8:
        H = .5
    else:
        with open(path['model_%s'%len(hydro_vector)], 'rb') as fin:
            model = pickle.load(fin)
        H = round(model.predict_proba(hydro_vector)[0][1], 3)
    return H


def GetRankExp(rank):
    return 1 / (1 + exp(5 * (rank - 2)))


def FxScore(Rm, Rn, H, RA, mismatch, comb):
    score = Rm * (1 - Rn / 2 ** mismatch) * RA * comb * H
    return score


def OneEpitScore(epit):
    allele, mtpep, mtpep_score, mtpep_aff, mtpep_rank, junc = epit.tolist()
    i = 0
    while i<3:
        i += 1
        try:
            wtpep, mismatch = GetWildpep(mut_pep=mtpep)
            wtpep_score, wtpep_aff, wtpep_rank = GetMHCpan(allele=allele, pep=wtpep)
            mtpep_tap, mtpep_cleavage, mtpep_comb = GetCTLpan(allele=allele, pep=mtpep)
            wtpep_tap, wtpep_cleavage, wtpep_comb = GetCTLpan(allele=allele, pep=wtpep)
            R = GetR(pep=mtpep)
            H = GetH(pep=mtpep)
            Rm = GetRankExp(float(mtpep_rank))
            Rn = GetRankExp(float(wtpep_rank))
            score = FxScore(Rm=Rm, Rn=Rn, H=float(H), RA=float(R), mismatch=float(mismatch), comb=float(mtpep_comb))
            if score < 0:
                score = 0
            out = '\t'.join(map(str, (allele, str(mismatch), mtpep, mtpep_score, mtpep_aff, mtpep_rank, mtpep_comb,
                         wtpep, wtpep_score, wtpep_aff, wtpep_rank, wtpep_comb, H, R, score, junc)))
            return out
        except Exception as e:
            logging.warning("Some error when process %s score, retry! error: %s" % (mtpep, e))


def ProcessNeo(process):
    logging.info("Begin Neoantign Score Process ...")
    with open(path['neo'], 'w') as fout:
        fout.write('HLA\tmismatch\t')
        fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
        fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
        fout.write('Hydro_Model\tR\tScore\tJunc_info\n')

        if os.path.getsize(path['epit']):
            epits = pd.read_table(path['epit'], header=None, sep='\s+')
            logging.info("Epit number is: %s" % len(epits))
            outs = MultiIters(its=epits, func=OneEpitScore, process=process)
            for out in outs:
                fout.write(out+'\n')

    neo = pd.read_table(path['neo'], header=0, sep='\s+')
    logging.info("Score info: sum %s, min %s, max %s, mean %s, median %s " % (
            neo.Score.sum(), neo.Score.min(), neo.Score.max(), neo.Score.mean(), neo.Score.median()))
    neo.sort_values(by='Score', ascending=False, inplace=True)
    neo.to_csv(path['neo_sorted'], sep='\t', header=True, index=False)


def CheckParameter(args):
    # check the parameter
    if not os.path.exists(args.junc):
        logging.error('Junction file "%s" not exists, please check!', args.junc)
        sys.exit(1)
    if not os.path.exists(args.genome):
        logging.error('Reference genome file "%s" not exists, please check!', args.genome)
        sys.exit(1)

    peplens = [int(l.strip()) for l in str(args.length).strip().split(',')]
    if not (7 < min(peplens) < 12 and 7 < max(peplens) < 12):
        logging.error('Peptide length set wrong, please check!')

    if args.bam:
        if not os.path.exists(args.bam):
            logging.error('Bam file "%s" not exists, please check!', args.bam)
            sys.exit(1)
        else:
            bai = os.path.join(os.path.dirname(os.path.abspath(args.bam)), os.path.basename(args.bam) + '.bai')
            if not os.path.exists(bai):
                logging.warning('Bam index file "%s" not exist, please check!', bai)
                sys.exit(1)

    with open(args.junc) as f:
        lines = f.readline().strip().split()
    if args.columns:
        columns = ['chrom', 'intron_start', 'intron_stop', 'unique_junction_reads']
        if not set(columns).issubset(set(lines)):
            logging.error('The junction file not contain necessary columns, please check!')
            sys.exit(1)
    else:
        if len(lines) != 9:
            logging.error('Junction file format error, please check if should set --columns')
            sys.exit(1)

    if not os.path.isdir(args.outdir):
        logging.warning('Out directory "%s" not exist, create!', args.outdir)
        os.makedirs(args.outdir)


def DefinePath(args):
    global path
    path = { }
    path['sdir'] = os.path.dirname(os.path.abspath(__file__)) # script dir

    path['bedtools'] = os.path.join(path['sdir'], 'src/bedtools')
    path['pepmatch'] = os.path.join(path['sdir'], 'src/pepmatch_db_x86_64')
    path['netMHCpan'] = os.path.join(path['sdir'], 'src/netMHCpan-4.0/netMHCpan')
    path['netCTLpan'] = os.path.join(path['sdir'], 'src/netCTLpan-1.1/netCTLpan')

    path['gtex'] = os.path.join(path['sdir'], 'data/Norm_SJ.tab')
    path['refseq_bed'] = os.path.join(path['sdir'], 'data/iRefSeq.bed')
    path['refseq_fa'] = os.path.join(path['sdir'], 'data/Norm_Protein.fasta')
    path['iedb'] = os.path.join(path['sdir'], 'data/iedb.fasta')

    for i in range(9,12):
        path['model_%s' % i] = os.path.join(path['sdir'], 'data/hg_xgb_%s.dat' % i)

    path['tdir'] = os.path.join(args.outdir, 'tmp') # tmp dir

    path['iso_bed'] = os.path.join(path['tdir'], 'putative_isoform.bed')
    path['iso_fa'] = os.path.join(path['tdir'], 'putative_isoform.fa')
    path['prot_fa'] = os.path.join(path['tdir'], 'putative_protein.fa')
    path['pep_txt'] = os.path.join(path['tdir'], 'putative_peptide.txt')
    path['affit'] = os.path.join(path['tdir'], 'netMHCpan.affinit')
    path['epit'] = os.path.join(path['tdir'], 'putative_epit.txt')
    path['neo'] = os.path.join(path['tdir'], 'putative_neo.txt')
    for i in range(8,12):
        path['refpep_%s'%i] = os.path.join(path['tdir'], 'reference_peptide_%s.txt' % i)

    path['neo_sorted'] = os.path.join(args.outdir, 'putative_neo.txt')


def main():
    print('\n' + "#" * 80)
    args = InputParser()
    logging.info('Start the program with %s\n', args)
    CheckParameter(args)
    DefinePath(args)
    if not os.path.exists(path['tdir']):
        os.mkdir(path['tdir'])

    ProcessSJ(args.junc, args.columns, args.reads, args.psi, args.bam, args.rpkm, args.process)
    ProcessIsoform(args.genome, args.length, args.allele, args.bind)
    ProcessNeo(args.process)
    shutil.rmtree(path['tdir'])
    logging.info('All Process Done!')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('keyboard interrupt')


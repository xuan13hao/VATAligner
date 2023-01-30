# pip install gtfparse, pandas, pysam, Biopython
# python3 gtf_extract.py match_tab SIRV_E2C.gtf 
from collections import defaultdict
from distutils.command.build import build
from operator import index
from os import lseek
from gtfparse import read_gtf
import numpy as np
import pandas as pd
import sys
import pysam
import re
from Bio import SeqIO


# >2:10252964-10252976
def readReadName(input_file): # truth
    truth = {}
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        seq_id = fasta.id
        seq_name = fasta.id
        m = re.match(r'(.*):(\d+)-(\d+)', seq_name)
        seq_name = m.group(1)
        s = int(m.group(2))
        e = int(m.group(3))
        # print(s,":",e)
        truth[seq_id] = (s,e) 
    # print(truth)
    print("Total reads num = ", len(truth))
    return truth

def buildSamPD(infile):
    samfile = pysam.AlignmentFile(infile,"r")
    colums = ['qseqid','qstart','qend','sseqid','sstart','send','Map_Qual'] 
    df_sam = []
    for read in samfile.fetch():
        if read.reference_name == None:
            qs = 0
            qe = 0
            rs = 0
            re = 0
        else:
            qs = read.query_alignment_start + 1
            qe = read.query_alignment_end + 1
            rs = read.reference_start + 1
            re = read.reference_end + 1
        record = [read.query_name,qs,qe,read.reference_name,rs,re,read.mapping_quality]
        df_sam.append(record)
    df = pd.DataFrame(df_sam, columns=colums)
    # print(len(df))
    t_dic = defaultdict(list)
    tabfile = df
    for i in tabfile.index:
        r = str(tabfile.loc[i,'sseqid'])
        q = tabfile.loc[i,'qseqid']
        # m = re.match(r'(.*):(\d+)-(\d+)', q)
        # r_seq_name = str(m.group(1))
        ss = tabfile.loc[i,'sstart']
        se = tabfile.loc[i,'send']
        # if r_seq_name == r:
        t_dic[q].append((int(ss),int(se)))
    # print(df)
    return t_dic

def readTabularFormatConvert(match_file):
    colums = ['qseqid','qstart','qend','sseqid','sstart','send','Map_Qual'] 
    df_tab = pd.read_csv(match_file,header=None,sep = '\t',names=colums)
    print(df_tab)
    t_dic = defaultdict(list)
    tabfile = df_tab
    for i in tabfile.index:
        r = str(tabfile.loc[i,'sseqid'])
        q = tabfile.loc[i,'qseqid']
        # m = re.match(r'(.*):(\d+)-(\d+)', q)
        # r_seq_name = str(m.group(1))
        ss = tabfile.loc[i,'sstart']
        se = tabfile.loc[i,'send']
        # if r_seq_name == r:
        t_dic[q].append((int(ss),int(se)))
    # print(t_dic)
    return t_dic


def readTabularFormat(match_file):
    columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    df_tab = pd.read_csv(match_file,header=None,sep = '\t',names=columns)
    t_dic = defaultdict(list)
    tabfile = df_tab
    for i in tabfile.index:
        r = str(tabfile.loc[i,'sseqid'])
        q = tabfile.loc[i,'qseqid']
        m = re.match(r'(.*):(\d+)-(\d+)', q)
        r_seq_name = str(m.group(1))
        ss = tabfile.loc[i,'sstart']
        se = tabfile.loc[i,'send']
        # if r_seq_name == r:
        t_dic[q].append((int(ss),int(se)))
    # print(t_dic)
    return t_dic


def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 > start2 and end2 > start1

def statistic(reads_truth, tabfile):
    num_exon = 0
    mapped_reads = 0
    mapped_count = {}
    mapped = {}
    for k in reads_truth:
        # print(k,":",v)
        if k in tabfile:
            (r1,r2) = reads_truth.get(k)
            algn_count = 0
            flag = 0
            for (p1,p2) in tabfile.get(k):
                if overlap(r1,r2,p1,p2):
                    algn_count = algn_count + 1
                    flag = 1
                    # print(p1,"-",p2)
                else:
                    flag = 0
                    algn_count = 0
            mapped_count[k] = algn_count
            mapped[k] = flag
        else:
            mapped_count[k] = 0
            mapped[k] = 0
        # else:
        #     unmapped[k] = unmapped[k] + 1
    for k, v in mapped.items():
        if v > 0:
            mapped_reads = mapped_reads + 1
    # for k, v in mapped_count.items():
    #     if v > 0:
    #         num_exon = num_exon + v
    print("mapped exon reads = ",mapped_reads)
    # print("exon count = ",num_exon)
    return num_exon

match_file = sys.argv[2]
fasta_file = sys.argv[1]
# print(readTabularFormat(match_file))
read_dict = readReadName(fasta_file)
# print(gtf_dict)
match = buildSamPD(match_file)
statistic(read_dict,match)
# for key in gdict:
#     print(key,':',gdict.get(key))
    # for a,b in val:
    #     print(a,b)
# print(grdth_dic)
        
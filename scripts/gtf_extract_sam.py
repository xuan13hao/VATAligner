# pip install gtfparse, pandas, pysam
# python3 gtf_extract.py alig_SIRV.sam SIRV_E2C.gtf 
from collections import defaultdict
from distutils.command.build import build
from gtfparse import read_gtf
import numpy as np
import pandas as pd
import sys
import pysam
import re
from Bio import SeqIO


# >ENSANIT00000000002.1_431_929_0:0:0_1:0:0_0/1
def readReadName(input_file): # truth
    reads_list = []
    colums = ['seqname','start','end']
    truth = defaultdict(list)
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        seq_name = fasta.id
        l = len(fasta.seq)
        m = re.match(r'(.*)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(.*)_(.*)_(.*)', seq_name)
        seq_name = m.group(1)
        s = int(m.group(2))
        e = int(m.group(2)) + l
        truth[seq_name].append((s,e)) 
    # df = pd.DataFrame(reads_list, columns=colums)
    return truth


def readGTFFormat(gtf_file):
    df_gtf = read_gtf(gtf_file)
    # print(df_gtf)
    gtf_dic = defaultdict(list)
    #signals = ['three_prime_utr','five_prime_utr','stop_codon','start_codon','Selenocysteine']
    for i in df_gtf.index:
        if df_gtf.loc[i,'feature'] == 'exon' and df_gtf.loc[i,'strand'] == '+':
            # t_id = df_gtf.loc[i,'transcript_id']+'.'+df_gtf.loc[i,'transcript_version']
            t_id = df_gtf.loc[i,'transcript_id']
            k = t_id
            ref_name = df_gtf.loc[i,'seqname']
            s = df_gtf.loc[i,'start']
            e = df_gtf.loc[i,'end']
            gtf_dic[k].append((int(s),int(e),ref_name))
    return gtf_dic

def getGroundTruth(gtf_dict,read_dict):
    grdth_dic = defaultdict(list)
    for key in read_dict:
        if key in gtf_dict:
            for p1,p2 in read_dict.get(key):
                # print("read locs: ",p1,"-",p2)
                cross_exon_dis = 0
                query_l = p1
                query_r = p2
                for p3,p4,ref in gtf_dict.get(key):
                    # print("grt locs: ",p3,"-",p4)
                    truth_l = p3
                    truth_r = p4
                    r1 = truth_l + query_l
                    r2 = truth_l + query_r
                    cross_exon_dis = r2 - truth_r
                    k = str(key)+"_"+str(ref)
                    if cross_exon_dis <= 0:
                        a = r1
                        b = r2
                        # print("non cross grt locs: ",a,"-",b)
                        grdth_dic[k].append((int(a),int(b)))
                        break
                    else:
                        a = r1
                        b = truth_r
                        # print("corss grt locs: ",a,"-",b)
                        grdth_dic[k].append((int(a),int(b)))
                        query_l = 0
                        query_r = cross_exon_dis
    num_grd = 0
    for key in grdth_dic:
        x = len(grdth_dic.get(key))
        num_grd = num_grd + x
    print("ground truth = ",num_grd)    
    return grdth_dic
        # print(key)
        # print(read_dict.get(key))
        # print(gtf_dict.get(key))

def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 > start2 and end2 > start1

def statistic(grd, tabfile):
    g_dic = defaultdict(list)
    t_dic = defaultdict(list)
    num_exon = 0
    exon_counts = []
    # for k,v in grd.items():
    #     print(k,":",v) 
    for i in tabfile.index:
        r = tabfile.loc[i,'sseqid']
        q = tabfile.loc[i,'qseqid']
        m = re.match(r'(.*)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(.*)_(.*)_(.*)', q)
        q_seq_name = m.group(1)
        # print(q_seq_name)
        ss = tabfile.loc[i,'sstart']
        se = tabfile.loc[i,'send']
        k = str(q_seq_name)+"_"+str(r)
        t_dic[k].append((int(ss),int(se)))
    for k,v in t_dic.items(): 
        print(k,":",v) 
    for key,val in grd.items():
        if key in t_dic:
            exon = []
            for (g1,g2) in val:
                exon_signal = 0
                for (t1,t2) in t_dic.get(key):
                    if overlap(g1,g2,t1,t2) and abs(t2-g1)/abs(g2-g1) > 0.8 or abs(g2-t1)/abs(g2-g1) > 0.8 or g1 < t1 or g2 > t2 :
                        exon_signal = 1
                exon.append(exon_signal)
            exon_counts.append(exon)
    for i in exon_counts:
        for j in i:
            num_exon = num_exon + j
    return num_exon


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
    # print(df)
    return df


match_file = sys.argv[3]
gtf_file = sys.argv[1]
fasta_file = sys.argv[2]
gtf_dict = readGTFFormat(gtf_file)
# print(readTabularFormat(match_file))
read_dict = readReadName(fasta_file)
# print(gtf_dict)
gdict = getGroundTruth(gtf_dict,read_dict)
match = buildSamPD(match_file)
exon = statistic(gdict,match)

print("mapped exon = " ,exon)
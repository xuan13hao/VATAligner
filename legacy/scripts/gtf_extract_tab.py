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
        m = re.match(r'(.*)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(.*)_(.*)_(.*)', seq_name)
        seq_name = m.group(1)
        #ENSATET00000014582.2
        # print(seq_name)
        # seq_name_ = re.match(r'(.*).(\d)', seq_name)
        # seq_name = seq_name_.group(1)
        l = len(fasta.seq)
        s = int(m.group(2))
        e = int(m.group(2)) + l
        # print(s,":",e)
        truth[seq_name].append((s,e)) 
        # truth[seq_name].append((int(m.group(2)),int(m.group(3)))) 
    # df = pd.DataFrame(reads_list, columns=colums)
    # print(df)
    print("Total reads num = ", len(truth))
    return truth



def readTabularFormat(match_file):
    columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    df_tab = pd.read_csv(match_file,header=None,sep = '\t',names=columns)
    # print("Total mapped num = ",len(df_tab))
    return df_tab

def readGTFFormat(gtf_file):
    df_gtf = read_gtf(gtf_file)
    # print(df_gtf['strand'])
    #print(df_gtf['seqname'])
    gtf_dic = defaultdict(list)
    #signals = ['three_prime_utr','five_prime_utr','stop_codon','start_codon','Selenocysteine']
    for i in df_gtf.index:
        #df_gtf.loc[i,'strand'] == '.'
        if df_gtf.loc[i,'feature'] == 'exon' and df_gtf.loc[i,'strand'] == '+':
            t_id = df_gtf.loc[i,'transcript_id']+'.'+df_gtf.loc[i,'transcript_version']
            # t_id = df_gtf.loc[i,'transcript_id']
            k = t_id
            ref_name = df_gtf.loc[i,'seqname']
            s = df_gtf.loc[i,'start']
            e = df_gtf.loc[i,'end']
            gtf_dic[k].append((int(s),int(e),ref_name))
    # print(gtf_dic)
    num_exon = 0
    for i in gtf_dic:
        num_exon = num_exon + len(gtf_dic.get(i))
    # print("Total GTF counts = ", len(gtf_dic))
    return gtf_dic

# chromosome:ASM359560v1:Chromosome:27247:28926:1 
# primary_assembly:fAnaTes1.2:10:18869194:18878717:-1 
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
                    # if r1 < p4 and r2 < p4:
                    #     a = r1
                    #     b = r2
                    #     grdth_dic[key].append((int(a),int(b),ref))
                    #     break
                    # if r1 < p4 and r2 > p4:
                    #     a = r1
                    #     b = p4
                    #     cross_exon_dis = p2 - p1 - (p4 - r1)
                    #     grdth_dic[key].append((int(a),int(b),ref))
                    # if cross_exon_dis > 0:
                    #     grdth_dic[key].append(int((p3),int(p3+cross_exon_dis),ref))  
    # print(grdth_dic)
    num_grd = 0
    for key in grdth_dic:
        x = len(grdth_dic.get(key))
        num_grd = num_grd + x
    print("grh = ",len(grdth_dic))
    # print("ground truth = ",num_grd)    
    return grdth_dic


def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 > start2 and end2 > start1

def statistic(grd, tabfile):
    result_dic = defaultdict(list)
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
        #ENSATET00000014582.2
        # seq_name_ = re.match(r'(.*).(\d)', q_seq_name)
        # q_seq_name = seq_name_.group(1)
        ss = tabfile.loc[i,'sstart']
        se = tabfile.loc[i,'send']
        k = str(q_seq_name)+"_"+str(r)
        t_dic[k].append((int(ss),int(se)))
    print("Mapped reads = ",len(t_dic))
    # for k,v in t_dic.items(): 
    #     print(k,":",v) 
    for key,val in grd.items():
        if key in t_dic:
            exon = []
            count_exon = 1
            for (g1,g2) in val:
                exon_signal = 0
                for (t1,t2) in t_dic.get(key):
                    #overlap(g1,g2,t1,t2) and abs(t2-g1)/abs(g2-g1) > 0.8 or abs(g2-t1)/abs(g2-g1) > 0.8 or g1 < t1 or g2 > t2
                    if overlap(g1,g2,t1,t2) :
                        exon_signal = 1
                exon.append(exon_signal)
            exon_counts.append(exon)
            result_dic[key] = count_exon
    for i in exon_counts:
        for j in i:
            num_exon = num_exon + j
    print("result len = ",len(result_dic))
    return num_exon

match_file = sys.argv[3]
gtf_file = sys.argv[1]
fasta_file = sys.argv[2]
gtf_dict = readGTFFormat(gtf_file)
# print(readTabularFormat(match_file))
read_dict = readReadName(fasta_file)
# print(gtf_dict)
gdict = getGroundTruth(gtf_dict,read_dict)
match = readTabularFormat(match_file)
exon = statistic(gdict,match)

print("mapped exon = " ,exon)
# for key in gdict:
#     print(key,':',gdict.get(key))
    # for a,b in val:
    #     print(a,b)
# print(grdth_dic)
        
import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'r')


# Usage
import argparse
from argparse import RawTextHelpFormatter
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
# m = re.match(r'(.*)_(\d+)_(\d+)_(.*)_(.*)_(.*)', seq_name)
# seq_name = m.group(1)

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Extract forward from a FASTA file generated by DWGSIM\n',
	usage='\n  %(prog)s [OPTIONS] [INPUTFASTA]  [OUTPUTFASTA]')
parser.add_argument('fasta', metavar='FASTA', help='input FASTA file', nargs='+')
parser.add_argument('--version', action='version', version='v0.2')
args = parser.parse_args()

# >CWI35_17925_256_245_1_0_0_0_0:0:0_2:0:0_0/1
# >(.*)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(.*)_(.*)_(.*)
sample = args.fasta[0]
outfile = args.fasta[1]
# fastqin = open(sample,"r")
# print(sample)
newSEQS = []
f1_signal = 150
for s in SeqIO.parse(sample, 'fasta'):
	newseq = s.seq
	len_ = len(newseq)
	if len_ > f1_signal:
		newseqREC = SeqRecord(newseq, id=s.id, description=s.description)
		newSEQS.append(newseqREC)

with open(outfile, 'w') as handle:
	for r in newSEQS:
		SeqIO.write(r, handle, 'fasta')
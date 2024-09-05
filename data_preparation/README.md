# Benchmark Data Preparation 
## Protein query and target datasets generation from Pfam for benchmark
The process to generate target.fasta and query.fasta from Pfam-A.fa involves two main steps. First, we randomly select a specific number of sequences from a large input FASTA file, such as the Pfam-A database, while imposing a limit on the number of sequences from the same Pfam family. This ensures diversity in the selected sequences by preventing any single family from dominating the dataset. The selected sequences are written to a new FASTA file, which serves as the basis for further analysis. In the second step, we split these selected sequences into two distinct files: a query set and a search database. A specified number of sequences are randomly chosen as query sequences and written to the query.fasta file. The remaining sequences are compiled into the target.fasta file, which will serve as the search database. 

## Short genomic datasets generation for benchmark
In the short-read genomic simulation study, Mason2 was employed to generate synthetic genomic reads extracted from a reference genome. Mason2 is a read simulator tool designed to replicate the sequencing process, allowing researchers to simulate short-read datasets with realistic sequencing errors, which are useful for benchmarking and testing various genomic tools. In this study, Mason2 was particularly useful for producing datasets with controlled characteristics, such as read length and sequencing error rates. The generated reads, formatted in the standard FASTQ format, are suitable for use as input in a wide range of short-read alignment tools, facilitating fair comparisons between them. For consistency, most of the synthetic datasets were extracted from the human genome (GRCh38), with Mason2 allowing for the generation of reads of different lengths and quantities as required for the study's objectives. By adjusting parameters such as the number of reads and the lengths of the sequences, Mason2 ensured that the simulated datasets were tailored to specific experimental conditions, https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_simulator.

```
mason2-2.0.9-Linux-x86_64/bin/mason_simulator -ir GCA_000001405.22_GRCh38.p7_genomic.fna -n 5000000 --illumina-read-length 150 -o hg38_reads.fa -oa hg38.reads.sam
```

## Long genomic datasets generation for benchmark
We simulated datasets from the human genome using Badread, a tool developed by Wick (2019), https://github.com/rrwick/Badread, specifically for generating synthetic reads that mimic the error profiles of real sequencing technologies. In our simulation, we used Badread to create reads with an average error rate of 12% and an average read length of 2,500 bases. The simulation was carried out using the cami_ref.fna reference genome. The command employed specified various parameters, including a random error model and ideal quality scores (--error_model random --qscore_model ideal), while avoiding the introduction of glitches, junk reads, random reads, and chimeras (--glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0). Additionally, we set the identity to vary around 30% with a standard deviation of 3%, and the read lengths were drawn from a normal distribution with a mean of 20,000 bases and a standard deviation of 1,000 bases (--identity 30,3 --length 20000,1000). No adapter sequences were added to the reads. The output of this simulation was compressed into a FASTQ file (cami_reads.fastq.gz) for further analysis, providing a robust dataset for benchmarking various bioinformatics tools and methodologies.
```
badread simulate --reference GRCh38_genomic.fn --quantity 1x --error_model random \
    --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    --identity 30,3 --length 20000,1000 --start_adapter_seq "" --end_adapter_seq "" \
    | gzip > GRCh38_genomic.long.reads.fastq.gz
```


## Short RNA-seq datasets generation for benchmark
In the process of short-reads simulation study , we applied a methodology inspired by Baruzzo et al. (2017) by http://bioinf.itmat.upenn.edu/BEERS/bp1/datasets.html, which emphasized simulation-based evaluation to assess the accuracy of aligners at various levels, including base, read, and exon junction accuracy. The method involved generating simulated RNA-seq data with varying levels of complexity and error rates using the BEERS, by https://github.com/itmat/beers_simulator/tree/kat, simulation engine. This approach allowed for a controlled assessment of how well different aligners handle challenges such as splice junctions, indels, and polymorphisms. By introducing high levels of complexity, the benchmark effectively highlighted the strengths and weaknesses of each aligner, particularly in scenarios involving difficult-to-align regions or cross-species data mapping. 

## Long RNA-seq datasets generation for benchmark

The simulation study was inspired by https://github.com/kkrizanovic/RNAseqEval, and make some difference and important to our study. The simulation workflow in our study was given here.

Synthetic datasets were created from following genomes and annotations:
  - Homo Sapiens GRCh38 (human) Reference genome, and annotations (version 94) from Ensemble.
 ---

We built two categories of simulated transcriptomes. The first category only has one simulated transcriptome (called as “H-all” transcriptome) which is from all the coding genes of human, and the second category has three simulated transcriptomes (called as “H-se”, “M-se” and “F-se” transcriptomes respectively) which are from three sets of randomly selected genes of human, mouse and fruit fly, respectively. Each of the transcriptomes was used to generate a series of simulated datasets with various sequencing models. The following are some details for the generation of “H-se”, “M-se” and “F-se” transcriptomes and datesets.

RNA-seq simulation datasets were generated by the following workflow:
  1. Download reference genome and annotation file with the same version information. The annotation file should contain the comprehensive gene annotations on the reference chromosomes only. This is the main annotation for most users. Or you can filter gene annotations on scaffolds, assembly patches and alternate loci if there exists.
  2. Extract genes with single splicing isoform, genes with alternative splicing isoforms and genes with short exons (< 30bp). For each kinds of genes, output all the corresponding annotations into three different files by suffix “_AS.gtf”, “_SS.gtf” and “_short.gtf” separately.
  3. Generate transcriptomes from pre-processed three kinds of annotations and reference genome. Filter transcripts short than 200 bps and combine generated transcriptomes together for the input of PBSIM and NanoSim.
  4. For  “PacBio ROI reads”, “PacBio subreads”, “ONT 2D reads”, “ONT 1D reads” datasets, simulate reads with different error rate (2%, 12%, 15%, 25%), and sequencing depth (4X, 10X, 30X) on the generated transcriptome by PBSIM. "PS-ONT reads” and “NS-ONT reads” were respectively generated by PBSim and NanoSim based on a real ONT dataset (SRA Accession Number: SRR2848544) to simulate more realistic ONT datasets.
 
### 1. Filter annotations not on chromosomes
As we know, the annotation file download from Ensemble contains all the gene annotations on chromosomes, scaffolds, assemble patches and alternate loci. But many of the scaffolds are unfinshed that should be removed. Or user can download file contains the comprehensive gene annotations on the reference chromosomes only which can be found at GENCODE.

### 2. Extract annotations of different gene types into groups
RNA-seq data is the products of gene transcription, but one gene may transcribe more than one isoforms due to alternative splicing of gene. What's more, an RNA aligner need to not only map exons to genome but also recognize isoforms of genes to reflect the gene structure. Thus, in order to reveal the performance of aligners mapping RNA-seq data to reference genome, we need to simulate RNA-seq data from different aspects to close to real data. In simulation, annotations for each species were seperated into three groups. The first group contains annotations for genes with single splicing isoform, the second group contains annotations for genes with alternative splicing isoforms, the third group contains annotations for genes with short exons (< 30bp). All annotations of each group will be output into files and suffixed by "_SS.gtf", "_AS.gtf" and "_short.gtf" separately. For alternative splicing genes, all isoforms will be extracted. Annotation grouping was down using Annotation_Load.py script in this repository, i.e.

```
python Annotation_grouping.py genome.gtf
```
### 3. Generating transcriptomes
Transcriptomes are generated from processed annotations and genome reference using the script generate_transcriptome.py (city from https://github.com/kkrizanovic/RNAseqEval) in this repository. Since the annotations were separated into three groups, a transcriptom (or a set of transcripts) was generated for each group.
```python
python generate_transcriptome.py SS.gtf ref.fa SS_transcriptome.fa
python generate_transcriptome.py AS.gtf ref.fa AS_transcriptome.fa
python generate_transcriptome.py short.gtf ref.fa short_transcriptome.fa
```
We combine the three transcriptome files together and then filter transcripts short than 200bp, this was done using fastqfilter.py script from https://github.com/isovic/samscripts with option minlen. The filtered transcriptome file is the input of PBSIM.
```
cat SS_transcriptome.fa AS_transcriptome.fa short_transcriptome.fa > merge_transcriptome.fa
python fastqfilter.py minlen 200 merge_transcriptome.fa > transcriptome_filter.fa
```

```
The PacBio subreads reads: 
	pbsim 	Transcripome_File (human) \
        --data-type CLR \
        --model_qc model_qc_clr	 \
        --length-mean 7800 \
        --length-min 100
        --difference-ratio 1:12:2 \
        --accuracy-mean 0.85 \
        --accuracy-min 0.75 \
        --depth 4/10/30
```
---
### One example for generation simulation dataset
Here is the pipeline of reads simulation, take human reference and annotations for example.

#### Get dataset

The following scripts give the pipeline of simulation data generation.
```
conda activate /home/h392x566/env/xuan
python Annotation_grouping.py Homo_sapiens.GRCh38.94.gtf
python generate_transcriptome.py Homo_sapiens.GRCh38.94_AS.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa GRch38.AS_transcriptome.fa
python generate_transcriptome.py Homo_sapiens.GRCh38.94_SS.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa GRch38.SS_transcriptome.fa
python generate_transcriptome.py Homo_sapiens.GRCh38.94_short.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa GRch38.short_transcriptome.fa

cat GRch38.AS_transcriptome.fa GRch38.SS_transcriptome.fa GRch38.short_transcriptome.fa >GRch38.merge_transcriptome.fa
pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim \
--data-type CLR --depth 40 \
--model_qc pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
--length-mean 3080 \
--length-sd 2211 \
--length-min 50 \
--length-max 50000 \
--accuracy-mean 0.95 \
--accuracy-sd 0.11 \
--accuracy-min 0.7 \
--difference-ratio 47:38:15 \
transcriptome_hg38.fa
```
## Chimeric datasets generation for benchmark
We simulated silico sequencing data for the benchmark experiment. We selected 18 microbial reference genomes  and used wgsim to generate two independent 50bp fragments and concatenated them together. We also inserted random bases (2-5bp) at the chimera junction to simulate potential crosslink artifacts and an error rate of 0.1% (Q30). 
```
wgsim -1 100 -2 100 -N mibcrobial.fna mibcrobial.r1.fq mibcrobial.r2.fq
cat mibcrobial.r1.fq mibcrobial.r2.fq> mibcrobial.fq
```
## Whole-Genome Alignment for benchmark

In our study, we used a benchmark method based on the MUMmer4 genome alignment system, as described by Marçais et al. (2018), to compare the human and chimpanzee genomes. This method was selected for its ability to handle large genomic datasets with high efficiency and precision, making it ideal for whole-genome alignments. The benchmark involved aligning the current assemblies of the human genome (GRCh38) and the chimpanzee genome (PanTro4), both of which are large and complex
### Dataset
```
ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/, human, GRCh38
www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001515.4/, chimpanzee
```
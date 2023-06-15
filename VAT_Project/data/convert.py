def convert_fastq_to_fasta(fastq_file1, fastq_file2, fasta_file):
    with open(fastq_file1, 'r') as f1, open(fastq_file2, 'r') as f2, open(fasta_file, 'w') as fw:
        while True:
            # Read four lines from each FASTQ file
            lines1 = [f1.readline().strip() for _ in range(4)]
            lines2 = [f2.readline().strip() for _ in range(4)]
            
            # Break the loop if end of file is reached in either file
            if not all(lines1) or not all(lines2):
                break
            
            # Extract the sequence identifier and sequence from the first file
            identifier = lines1[0][1:]
            sequence = lines1[1]
            
            # Write the sequence in FASTA format to the output file
            fw.write(f'>{identifier}\n{sequence}\n')
            
            # Skip the three lines in the second file (quality score, +, and quality scores)
            _ = [f2.readline() for _ in range(3)]

# Specify the file paths
fastq_file1 = 'r_1.fastq'
fastq_file2 = 'r_2.fastq'
fasta_file = 'combined.fa'

# Convert the FASTQ files to FASTA
convert_fastq_to_fasta(fastq_file1, fastq_file2, fasta_file)

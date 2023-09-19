import matplotlib.pyplot as plt
import sys
import pandas as pd

fasta_file = sys.argv[1]
# Read the BLAST tabular format file
blast_data = pd.read_csv(fasta_file, sep='\t', header=None)

# Select relevant columns: query start, query end, subject start, subject end
query_start = blast_data[6]
query_end = blast_data[7]
subject_start = blast_data[8]
subject_end = blast_data[9]

# Create a scatter plot with alignment diagonal lines
plt.figure(figsize=(8, 6))

for i in range(len(blast_data)):
    x = [query_start[i], query_end[i]]
    y = [subject_start[i], subject_end[i]]
    plt.plot(x, y, marker='o')

plt.xlabel("Query Position")
plt.ylabel("Subject Position")
plt.title("Alignment Diagonal in BLAST Hits")
plt.grid(True)
plt.show()

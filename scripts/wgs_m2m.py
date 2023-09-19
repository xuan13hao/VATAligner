
import sys

match_file = sys.argv[1]
blasttab_file = match_file

alignment_data = {}  # Dictionary to store alignment data

# Parse the BLAST tabular file
with open(blasttab_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        query_name = fields[0]
        subject_name = fields[1]
        query_start = int(fields[6])
        query_end = int(fields[7])
        subject_start = int(fields[8])
        subject_end = int(fields[9])
        alignment_identity = float(fields[2])
        alignment_score = float(fields[11])

        if query_name not in alignment_data:
            alignment_data[query_name] = []

        alignment_data[query_name].append((subject_name, query_start, query_end, subject_start, subject_end, alignment_identity, alignment_score))

average_identities = []

# Identify many-to-many best alignments and calculate average identity
for query_name, alignments in alignment_data.items():
    unique_alignments = set()
    non_overlapping_alignments = []
    
    # Sort alignments by identity in descending order
    alignments.sort(key=lambda x: x[5], reverse=True)

    for alignment in alignments:
        alignment_key = (alignment[1], alignment[2], alignment[3], alignment[4])
        if alignment_key not in unique_alignments:
            unique_alignments.add(alignment_key)
            non_overlapping_alignments.append(alignment)

    # Calculate average identity for many-to-many alignments
    total_identity = sum(alignment[5] for alignment in non_overlapping_alignments)
    average_identity = total_identity / len(non_overlapping_alignments)
    average_identities.append(average_identity)

# Calculate the overall average identity
overall_average_identity = sum(average_identities) / len(average_identities)

print("Overall Average Identity for Many-to-Many Alignments:", overall_average_identity)








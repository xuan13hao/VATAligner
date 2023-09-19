import sys
def parse_blasttab(file_path):
    alignments = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            query_id = parts[0]
            subject_id = parts[1]
            percent_identity = float(parts[2])
            alignment_length = int(parts[3])
            e_value = float(parts[10])
            
            if query_id not in alignments:
                alignments[query_id] = []
            
            alignments[query_id].append((subject_id, percent_identity, e_value))
    
    return alignments

def compute_best_alignments(alignments):
    best_alignments = {}
    
    for query_id, alignments_list in alignments.items():
        best_alignment = max(alignments_list, key=lambda x: x[1])  # Choose the alignment with highest percent identity
        best_alignments[query_id] = best_alignment
    
    return best_alignments

def calculate_average_identity(best_alignments):
    total_identity = sum(align[1] for align in best_alignments.values())
    average_identity = total_identity / len(best_alignments)
    return average_identity

# match_file = sys.argv[2]
fasta_file = sys.argv[1]
# Replace 'blast_results.txt' with your actual BLAST tabular file
blast_file_path = fasta_file

# Step 1: Parse the BLAST tabular file
alignments = parse_blasttab(blast_file_path)

# Step 2: Compute 1-to-1 best alignments
best_alignments = compute_best_alignments(alignments)

# Step 3: Calculate average identity
average_identity = calculate_average_identity(best_alignments)
print("Number :", len(alignments))
print("Best Alignment Number :", len(best_alignments))
print("Average Identity:", average_identity)

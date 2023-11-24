from Bio import SeqIO
import sys
def load_ground_truth(fasta_file):
    """
    Load the ground truth sequence IDs from a FASTA file and store them in a dictionary.
    The keys are sequence IDs, and the values can be any relevant information.
    """
    ground_truth = {}
    num = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        num = num + 1
        # Assuming you have additional information associated with the sequence IDs
        ground_truth[record.id] = record.description  # You can use any relevant data
        # print(record.id,":",record.description)
    return ground_truth,num

def compute_tp_fp_fn(blasttab_file, ground_truth):
    """
    Compute true positives (TP), false positives (FP), and false negatives (FN) by comparing
    sequence IDs in the BLASTTab file to the ground truth dictionary.
    """
    true_positives = 0
    false_positives = 0
    false_negatives = 0

    with open(blasttab_file, "r") as result_file:
        for line in result_file:
            fields = line.split("\t")
            # print(fields)
            if len(fields) >= 2:
                query_id = fields[0]
                reference_id = fields[1]
                # print(query_id, ":",reference_id)
                if query_id == reference_id and query_id in ground_truth:
                    true_positives += 1
                elif query_id != reference_id and (query_id in ground_truth or reference_id in ground_truth):
                    false_positives += 1
                else:
                    false_negatives += 1

    return true_positives, false_positives, false_negatives
# Define the paths to your ground truth FASTA file and the query results file (BLASTTab format).
ground_truth_file = sys.argv[1]
blasttab_file = sys.argv[2]

# Load the ground truth sequence IDs and associated information into a dictionary
ground_truth,num = load_ground_truth(ground_truth_file)
print(num)
# Compute TP, FP, and FN
true_positives, false_positives, false_negatives = compute_tp_fp_fn(blasttab_file, ground_truth)


# Print or use the results as needed
print("True Positives:", true_positives)
print("False Positives:", false_positives)
print("False Negatives:", false_negatives)


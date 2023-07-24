import sys
def read_fasta_file(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_sequence = ''
        for line in f:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences
def is_overlap(start1, end1, start2, end2):
    """
    Check if there is an overlap between two intervals.

    Parameters:
        start1 (int): The start position of the first interval.
        end1 (int): The end position of the first interval.
        start2 (int): The start position of the second interval.
        end2 (int): The end position of the second interval.

    Returns:
        bool: True if there is an overlap, False otherwise.
    """
    return not (end1 < start2 or start1 > end2)

def extract_start_end_positions(interval_str):
    """
    Extract the start and end positions from an interval string and convert them to integers.

    Parameters:
        interval_str (str): The interval string in the format "start-end".

    Returns:
        tuple: A tuple containing the start and end positions as integers.
    """
    interval_parts = interval_str.split('_')[0]
    start, end = interval_parts.split(':')[1].split('-')
    return int(start), int(end)


def count_mapped_reads(blasttab_file, simulated_reads):
    mapped_reads = set()
    with open(blasttab_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = fields[0]
            ref_s,ref_e = extract_start_end_positions(query_id)
            # print(query_id,"",fields[3])
            read_length  = int(fields[3])
            #is_overlap(ref_s,ref_e,int(fields[8]),int(fields[9]))
            if int(read_length) > 80 and is_overlap(ref_s,ref_e,int(fields[8]),int(fields[9])):
                mapped_reads.add(query_id)
    num_simulated_reads = len(simulated_reads)
    num_mapped_reads = len(mapped_reads)
    return num_mapped_reads, num_simulated_reads
#python3 spliced_benchmark.py Anabas_exons_forward_sim.fa Anabas_match_sim_fd_2.0 
if __name__ == "__main__":
    match_file = sys.argv[2]
    fasta_file = sys.argv[1]
    simulated_reads_fasta = fasta_file
    blasttab_file = match_file

    simulated_reads = read_fasta_file(simulated_reads_fasta)
    mapped_read_count, total_simulated_reads = count_mapped_reads(blasttab_file, simulated_reads)

    print(f"Total simulated reads: {total_simulated_reads}")
    print(f"Number of reads mapped to reference genome: {mapped_read_count}")
    print(f"Percentage of mapped reads: {(mapped_read_count / total_simulated_reads) * 100:.2f}%")

import sys
def partition_fasta(input_file, num_chunks):
    sequences = {}
    current_sequence_name = None
    current_sequence = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if current_sequence_name:
                    sequences[current_sequence_name] = current_sequence

                current_sequence_name = line[1:]
                current_sequence = [line]
            else:
                current_sequence.append(line)

    if current_sequence_name:
        sequences[current_sequence_name] = current_sequence

    chunk_size = len(sequences) // num_chunks
    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = (chunk_idx + 1) * chunk_size

        chunk_sequences = list(sequences.keys())[start_idx:end_idx]

        output_file = f"chunk_{chunk_idx + 1}.fasta"
        with open(output_file, 'w') as out:
            for seq_name in chunk_sequences:
                out.write('\n'.join(sequences[seq_name]) + '\n')

if __name__ == '__main__':
    fasta_file = sys.argv[1]  # Replace with your input FASTA file
    num_chunks = 100

    partition_fasta(fasta_file, num_chunks)



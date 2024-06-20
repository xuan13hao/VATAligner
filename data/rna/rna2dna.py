import sys
def rna_to_dna(rna_sequence):
    # RNA to DNA conversion: replace 'U' and 'u' with 'T'
    return rna_sequence.replace('U', 'T').replace('u', 'T')

def rna_to_dna_sequences(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            # Ensure the line is stripped of leading/trailing whitespace
            stripped_line = line.strip()
            # Process each line that is not a header (starting with '>')
            if not stripped_line.startswith('>'):
                dna_sequence = rna_to_dna(stripped_line)
                outfile.write(dna_sequence + '\n')
            else:
                # Write the header line unchanged
                outfile.write(line)

# Example usage
input_filename = sys.argv[1]
output_filename = sys.argv[2]

# Convert RNA sequences to DNA sequences and save to file
rna_to_dna_sequences(input_filename, output_filename)

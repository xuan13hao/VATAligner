sequence = "AAATTCCAAGCCCGATTGATAATTTCCATAATACAAGAAATCGTCTGATGAAGTACACAAGTCATCATTGTAGTAGCATGTAGTGGCGGTGGTACACTAA"  # Example sequence

complementary_mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}

complementary_sequence = "".join(complementary_mapping[base] for base in sequence)
print("Original sequence:", sequence)
print("Complementary sequence:", complementary_sequence)


def complementary_reverse_sequence(sequence):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reversed_sequence = sequence[::-1]
    complementary_sequence = ''.join(complement[base] for base in reversed_sequence)
    return complementary_sequence

print(complementary_reverse_sequence(sequence))
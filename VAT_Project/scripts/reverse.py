sequence = "TTCTTCATGCTGCTCCTCTGATCGCCGGTGAAAGTTGGTTTCATTGTAAT"  # Example sequence

complementary_mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}

complementary_sequence = "".join(complementary_mapping[base] for base in sequence)
print("Original sequence:", sequence)
print("Complementary sequence:", complementary_sequence)
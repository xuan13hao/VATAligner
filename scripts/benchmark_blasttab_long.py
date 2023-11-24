import sys
import re
def read_fasta(fasta_file):
    sequences = {}
    current_sequence = None

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    new_dict = {extract_uuid(new_key): extract_positions(new_key) for new_key, value in sequences.items()}
    return new_dict
def extract_positions(input_string):
    # Define a regular expression pattern to match the desired substring
    pattern = r',(\d+-\d+)'

    # Use re.search to find the first match of the pattern in the input string
    match = re.search(pattern, input_string)

    if match:
        # Extract the matched substring
        matched_substring = match.group(1)
        return matched_substring
    else:
        # If no match is found, return None or handle the case accordingly
        return None

def extract_uuid(input_string):
    # Split the input string by spaces
    parts = input_string.split()
    
    # Check if there are enough parts in the split string
    if len(parts) >= 1:
        # The first part should be the UUID
        uuid = parts[0]
        return uuid
    
    # If the input string doesn't have the expected format, return None or raise an exception
    return None
def extract_numbers(input_string):
    # Split the input string using the hyphen as a delimiter
    numbers = input_string.split('-')

    if len(numbers) == 2:
        # Convert the extracted substrings to integers
        try:
            num1 = int(numbers[0])
            num2 = int(numbers[1])
            return num1, num2
        except ValueError:
            return None, None
    else:
        # Return None if the input string doesn't contain two numbers separated by a hyphen
        return None, None
def read_blasttab(blasttab_file):
    mappings = []

    with open(blasttab_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            query_id, subject_id, percent_identity, alignment_length, mismatches, gap_openings, q_start, q_end, s_start, s_end, e_value, bit_score = line
            mappings.append({
                'query_id': query_id,
                'subject_id': subject_id,
                'q_start': int(q_start),
                'q_end': int(q_end),
                's_start': int(s_start),
                's_end': int(s_end),
            })

    return mappings
def is_overlap(start1, end1, start2, end2, min_overlap_percentage=80):
    """
    Check if there is an overlap between the first interval and the second interval with a specified minimum overlap percentage.

    Parameters:
        start1 (int): The start position of the first interval.
        end1 (int): The end position of the first interval.
        start2 (int): The start position of the second interval.
        end2 (int): The end position of the second interval.
        min_overlap_percentage (int): The minimum required overlap percentage (default is 80).

    Returns:
        bool: True if the overlap percentage is greater than or equal to min_overlap_percentage, False otherwise.
    """
    # Calculate the overlap length (intersection)
    overlap = max(0, min(end1, end2) - max(start1, start2))

    # Calculate the length of the first interval
    length1 = end1 - start1

    # Calculate the overlap percentage for the first interval
    overlap_percentage = (overlap / length1) * 100

    # Check if the overlap percentage is greater than or equal to min_overlap_percentage
    return overlap_percentage >= min_overlap_percentage


def calculate_accuracy(ground_truth, mappings):
    correct_mappings = {}
    total_mappings = len(ground_truth)

    for mapping in mappings:
        query_id = mapping['query_id']
        s_start = mapping['s_start']
        s_end = mapping['s_end']
        # ground_truth_seq = ground_truth[query_id]
        # mapped_seq = ground_truth_seq[s_start - 1:s_end]
        # Check if the string is a key in the dictionary
        if query_id in ground_truth:
            r_s, r_e = extract_numbers(ground_truth[query_id])
            if is_overlap(r_s,r_e,s_start,s_end):
                correct_mappings[query_id] = 1
    accuracy = (len(correct_mappings) / total_mappings) * 100
    # print(correct_mappings)
    return accuracy

if __name__ == "__main__":
    blasttab_file = sys.argv[2]
    fasta_file = sys.argv[1]
    ground_truth = read_fasta(fasta_file)
    # print(ground_truth)
    mappings = read_blasttab(blasttab_file)
    accuracy = calculate_accuracy(ground_truth, mappings)

    print(f"Accuracy: {accuracy:.2f}%")

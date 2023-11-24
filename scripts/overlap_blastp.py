import sys
def parse_blasttab_file(filename):
    """
    Parse a BLASTTAB format file and return a dictionary with (query_id, subject_id, subject_start, subject_end)
    as keys and a list of tuples containing (query_start, query_end, subject_start, subject_end) as values.
    """
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query_id = fields[0]
                subject_id = fields[1]
                query_start = int(fields[6])
                query_end = int(fields[7])
                subject_start = int(fields[8])
                subject_end = int(fields[9])

                key = (query_id, subject_id, subject_start, subject_end)

                if key not in data:
                    data[key] = []
                data[key].append((query_start, query_end, subject_start, subject_end))

    return data

def calculate_overlap_rate(file1, file2):
    """
    Calculate the overlap rate between two BLASTTAB format files based on the specified conditions
    using the dictionary created from file1.
    """
    data1 = parse_blasttab_file(file1)
    data2 = parse_blasttab_file(file2)

    overlap_count = 0
    total_count = 0

    for key in data2:
        if key in data1:
            for entry1 in data1[key]:
                for entry2 in data2[key]:
                    query_overlap = max(0, min(entry1[1], entry2[1]) - max(entry1[0], entry2[0]) + 1)
                    subject_overlap = max(0, min(entry1[3], entry2[3]) - max(entry1[2], entry2[2]) + 1)

                    if query_overlap > 0 and subject_overlap > 0:
                        overlap_count += 1

                total_count += 1

    overlap_rate = overlap_count / total_count if total_count > 0 else 0
    return overlap_count

# Example usage:
ground_truth_file = sys.argv[1]
blasttab_file = sys.argv[2]
overlap_count = calculate_overlap_rate(ground_truth_file, blasttab_file)
print(overlap_count)


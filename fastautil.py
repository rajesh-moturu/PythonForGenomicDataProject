"""
Author: Rajesh Moturu
Process fasta file into id, seq dictionary and get basic information like
#records (or sequences)
longest and shortest sequences along with its id
"""
def read_fasta(fasta):
    """
    Read fasta into dictionary as key value pairs
    """
    sequences = {}
    id, sequence = "", ""
    for line in fasta:
        line = line.rstrip()
        if line.startswith('>'):
            if id == "":
                header = line.split('|')[-1]
                id = header.split(' ')[0]
                continue
            sequences[id] = sequence
            header = line.split('|')[-1]
            id = header.split(' ')[0]
            sequence = ""
        else:
            sequence += line
    sequences[id] = sequence

    return sequences


def count_records(sequences):
    """
    count number of records a sequences dictionary has
    """
    return len(sequences.keys())


def find_longest_and_shortest_seq(sequences):
    """
    find longest and shortest sequences of a given sequences dictionary
    """
    seq_len = {}
    for id, seq in sequences.items():
        seq_len[id] = len(seq)

    sorted_seq = dict(sorted(seq_len.items(), key=lambda item: item[1]))

    max_key, min_key = list(sorted_seq.keys())[-1], list(sorted_seq.keys())[0]
    max_value, min_value = list(sorted_seq.values())[-1], list(sorted_seq.values())[0]

    max_seq, min_seq = {}, {}
    for key, value in sorted_seq.items():
        if value == max_value:
            max_seq[key] = value
        if value == min_value:
            min_seq[key] = value

    return max_seq, min_seq

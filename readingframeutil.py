"""
Author: Rajesh Moturu

In molecular biology, a reading frame is a way of dividing the DNA sequence of
nucleotides into a set of consecutive, non-overlapping triplets (or codons).
Depending on where we start, there are six possible reading frames: three in the
forward (5' to 3') direction and three in the reverse (3' to 5'). For instance,
the three possible forward reading frames for the sequence
AGGTGACACCGCAAGCCTTATATTAGC are:

AGG TGA CAC CGC AAG CCT TAT ATT AGC
A GGT GAC ACC GCA AGC CTT ATA TTA GC
AG GTG ACA CCG CAA GCC TTA TAT TAG C

These are called reading frames 1, 2, and 3 respectively. An open reading frame (ORF)
is the part of a reading frame that has the potential to encode a protein. It starts
with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA). For instance,
ATGAAATAG is an ORF of length 9.
"""

def get_all_orfs(sequences, id = None):
    '''
    Takes the id, sequence of all the sequences of a fasta file and returns ORFs
    in all the sequences along with its pos of occurence and id.
    if id is provided, only returns ORFs for that specific id.
    '''
    START_CODON = "ATG"
    STOP_CODONS = ["TAA", "TAG", "TGA"]

    # reading frame 1, 2, or 4
    def get_pos_and_orfs_per_seq(reading_frame, start_pos):
        orf = ""
        all_orfs = []
        start_codon_found, stop_codon_found = False, False
        for idx in range(start_pos, len(reading_frame), 3):
            codon = reading_frame[idx:idx+3]
            if codon == START_CODON:
                start_codon_found = True
                start_codon_pos = idx + 1  #  to know where the longest orf occured
            if codon in STOP_CODONS: stop_codon_found = True

            if not start_codon_found:
                continue

            orf += codon
            if stop_codon_found:
                all_orfs.append((start_codon_pos, orf))
                start_codon_found, stop_codon_found, orf = False, False, ""

        return all_orfs

    def get_orfs_from_all_sequences(start_pos):
        orfs_all_sequences = {}
        for _id, _seq in sequences.items():
            # if we need operation only for one ID, skip others
            if id is not None and _id != id:
                continue
            pos_and_orfs_per_seq = get_pos_and_orfs_per_seq(_seq, start_pos)
            if len(pos_and_orfs_per_seq) != 0:
                orfs_all_sequences[_id] = pos_and_orfs_per_seq

        return orfs_all_sequences

    return get_orfs_from_all_sequences


def get_longest_orf_per_each_seq(all_orfs, fall_back = ""):
    """
    Takes a dictionary of type:
    {id: [(pos, orf), (pos2, orf2), ...], ...}
    and returns the id and longest orf
    """
    id_wise_longest_orf = {}
    for id, pos_orf_list in all_orfs.items():
        orfs = [(pos, orf) for pos, orf in pos_orf_list]
        longest_orf = max(orfs, key=lambda item: len(item[1]))
        id_wise_longest_orf[id] = longest_orf

    return id_wise_longest_orf

def get_longest_orf_of_all_seq(id_wise_longest_orf):
    """
    Takes a dictionary of type {id: (pos, longest_orf), ...}
    and returns the id, pos and longest_orf of all
    """
    max_length = 0
    max_orf_id = ""
    for id, pos_orf_value in id_wise_longest_orf.items():
        pos, orf = pos_orf_value
        if max_length < len(orf):
            max_orf_id = id

    return id, id_wise_longest_orf[id]

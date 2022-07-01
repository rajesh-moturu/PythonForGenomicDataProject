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
    def get_pos_and_orfs_per_seq(sequence, start_pos):
        orfs = []
        for i in range(start_pos, len(sequence), 3):
            codon = sequence[i:i+3]
            # once we hit start codon, we will loop over until we hit first stop codon
            if codon == START_CODON:
                for j in range(i+3, len(sequence), 3):
                    codon = sequence[j:j+3]
                    if codon in STOP_CODONS:
                        orfs.append((i+1, sequence[i:j+3]))  # i+1 since pos is not index
                        break
        return orfs

    def get_orfs_from_all_sequences(reading_frame):
        orfs_all_sequences = {}
        for _id, _seq in sequences.items():
            # if we need operation only for one ID, skip others
            if id is not None and _id != id:
                continue
            pos_and_orfs_per_seq = get_pos_and_orfs_per_seq(_seq, reading_frame - 1)
            if len(pos_and_orfs_per_seq) != 0:
                orfs_all_sequences[_id] = pos_and_orfs_per_seq

        return orfs_all_sequences

    return get_orfs_from_all_sequences


def get_longest_orf_per_each_seq(all_orfs):
    """
    Takes a dictionary of type:
    {id: [(pos, orf), (pos2, orf2), ...], ...}
    and returns the id and longest orf
    """
    id_wise_longest_orf = {}
    for id, pos_orf_list in all_orfs.items():
        longest_orf_pos, longest_orf = pos_orf_list[0]
        for pos, orf in pos_orf_list[1:]:
            if len(orf) > len(longest_orf):
                longest_orf_pos, longest_orf = pos, orf
        id_wise_longest_orf[id] = (longest_orf_pos, longest_orf)

    return id_wise_longest_orf

def get_longest_orf_of_all_seq(id_wise_longest_orf):
    """
    Takes a dictionary of type {id: (pos, longest_orf), ...}
    and returns the id, pos and longest_orf of all
    """
    lg_orf_id = list(id_wise_longest_orf.keys())[0]
    lg_orf_pos, lg_orf = id_wise_longest_orf[lg_orf_id]
    for id, (pos, orf) in id_wise_longest_orf.items():
        if len(orf) > len(lg_orf):
            lg_orf_id, lg_orf_pos, lg_orf = id, pos, orf

    return lg_orf_id, (lg_orf_pos, lg_orf)

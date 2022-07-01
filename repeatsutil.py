"""
Author: Rajesh M

Basic operations on repeats of a given dictionary of id, sequences such
#count of occurence of each repeat
#repeats of a given length 'n'

A repeat is a substring of a DNA sequence that occurs in multiple copies (more
than one) somewhere in the sequence. Although repeats can occur on both the forward
and reverse strands of the DNA sequence, we will only consider repeats on the forward
strand here.
"""
from collections import Counter

def get_all_repeats(sequences, n):

    def count_repeats_per_seq(sequence):
        n_mers = [sequence[i:i+n] for i in range(len(sequence))]
        counts = dict(Counter(n_mers))
        repeats = dict([(n_mer, count) for n_mer, count in counts.items() if count > 1])
        return repeats

    def make_dict_with_repeats():
        repeats = {}
        for id, seq in sequences.items():
            repeats_per_seq = count_repeats_per_seq(seq)
            # add the repeats occ to get most occ repeat in all seq.
            for repeat, counter in repeats_per_seq.items():
                repeats[repeat] = repeats[repeat] + counter if repeat in repeats.keys() else counter
        return repeats

    def most_freq_repeat_all_seq():
        repeats = make_dict_with_repeats()
        # print(sorted(repeats.items()))
        most_freq_repeat, most_count = "", 0
        check = ["CGCGCCG", "CATCGCC", "TGCGCGC", "AATGGCA"]
        for n_mer, count in repeats.items():
            if n_mer in check:
                print(f"repeat: {n_mer}, count: {count}")
            if count > most_count:
                most_freq_n_mer, most_count = n_mer, count
        return most_freq_n_mer, most_count

    return most_freq_repeat_all_seq

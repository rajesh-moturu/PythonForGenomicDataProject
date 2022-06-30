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

    def get_repeats(sequence):
        n_mers = [seq[i:i+n] for i in range(len(seq))]
        counts = dict(Counter(n_mers))
        repeats = dict([(n_mer, count) for n_mer, count in counts.items() if count > 1])
        return repeats

    for id, seq in sequences.items():
        repeats = get_repeats(seq)
        sorted_repeats = dict(sorted(repeats.items(), key=lambda item: item[1]))
        # get max repeat
        repeat_pairs = sorted_repeats.items()
        repeat_iter = iter(repeat_pairs)
        yield next(repeat_iter)

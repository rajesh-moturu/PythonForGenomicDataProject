import sys
import fastautil as fu
import readingframeutil as rf
import repeatsutil as rp

def main():
    if len(sys.argv) != 2:
        print("Need input fasta file.")
        sys.exit(1)
    fasta = open(sys.argv[1])

    print("\n======== RECORD COUNT DETAILS =======\n")
    # read counts
    sequences = fu.read_fasta(fasta)
    total_records = fu.count_records(sequences)
    print("Records count: ", total_records)

    # longest and shortest sequences
    print("\n======== LONGEST AND SHORTEST SEQ DETAILS =======\n")
    max_seqs, min_seqs = fu.find_longest_and_shortest_seq(sequences)
    print(f"Max sequences: {max_seqs}")
    print(f"Min sequences: {min_seqs}")

    # ORFs
    reading_frame = int(input("Enter the reading frame: "))
    get_orfs_from_all_sequences = rf.get_all_orfs(sequences, id="gi|142022655|gb|EQ086233.1|16")
    all_orfs = get_orfs_from_all_sequences(reading_frame)
    seq_wise_longest_orf = rf.get_longest_orf_per_each_seq(all_orfs)

    print("\n======== LONGEST ORF DETAILS =======\n")
    longest_orf_id, (start_pos, orf) = rf.get_longest_orf_of_all_seq(seq_wise_longest_orf)
    print(f"id: {longest_orf_id}")
    print(f"Pos: {start_pos}")
    print("Longest orf len: ", len(orf))

    print("\n======== REPEATS INFO =======\n")
    n = int(input("Please enter repeat length(n):\n"))
    most_freq_repeat_all_seq = rp.get_all_repeats(sequences, n)
    print(f"Most freq repeats:\n{most_freq_repeat_all_seq()}")

if __name__ == "__main__":
    main()

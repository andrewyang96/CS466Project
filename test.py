from algorithms.global_alignment import needleman_wunsch
from algorithms.local_alignment import smith_waterman

from Bio import SeqIO
import sys
sys.setrecursionlimit(25000)

if __name__ == '__main__':
    seqs = SeqIO.parse('dna.fasta', 'fasta')
    seq1 = next(seqs)[:1000]
    seq2 = next(seqs)[:1000]
    seq3 = next(seqs)[:100]

    global_result = needleman_wunsch(seq1, seq2)
    local_result = smith_waterman(seq1, seq3)

"""Global alignment using the Needleman–Wunsch algorithm.

Source: http://schulte.faculty.unlv.edu//BIO480/SequenceAlign.pdf
"""

import numpy as np

NONE, LEFT, UP, DIAG = range(4)

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap_start_penalty=0, gap_extend_penalty=-2):
    """Returns the two alignments and the optimal score."""
    assert match >= 0
    assert mismatch <= 0
    assert gap_start_penalty <= 0
    assert gap_extend_penalty <= 0

    m, n = len(seq1), len(seq2)
    scores = np.zeros((m+1, n+1))
    gap_i = np.ones((m+1,), dtype=bool)
    gap_j = np.ones((n+1,), dtype=bool)
    pointer = np.zeros((m+1, n+1))

    pointer[0, 0] = NONE
    scores[0, 0] = 0
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP
    scores[0, 1:] = gap_start_penalty + gap_extend_penalty * np.arange(n)
    scores[1:, 0] = gap_start_penalty + gap_extend_penalty * np.arange(m)

    gap_i[0] = False
    for i in range(1, m+1):
        gap_j[0] = False
        for j in range(1, n+1):
            gap_j[j] = True
            diag_score = scores[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up_score = scores[i-1, j] + (gap_extend_penalty if gap_i[i-1] else gap_start_penalty)
            left_score = scores[i, j-1] + (gap_extend_penalty if gap_j[j-1] else gap_extend_penalty)
            if diag_score >= up_score:
                if diag_score >= left_score:
                    scores[i, j] = diag_score
                    pointer[i, j] = DIAG
                    gap_i[i] = False
                    gap_j[j] = False
                else:
                    scores[i, j] = left_score
                    pointer[i, j] = LEFT
            else:
                if up_score >= left_score:
                    scores[i, j] = up_score
                    pointer[i, j] = UP
                else:
                    scores[i, j] = left_score
                    pointer[i, j] = LEFT

    aligned1 = []
    aligned2 = []
    while pointer[i, j] != NONE:
        if pointer[i, j] == DIAG:
            i -= 1
            j -= 1
            aligned1.append(seq1[i])
            aligned2.append(seq2[j])
        elif pointer[i, j] == LEFT:
            j -= 1
            aligned1.append('-')
            aligned2.append(seq2[j])
        elif pointer[i, j] == UP:
            i -= 1
            aligned1.append(seq1[i])
            aligned2.append('-')
        else:
            raise RuntimeError('NONE unexpectedly encountered')
    return ''.join(aligned1[::-1]), ''.join(aligned2[::-1]), scores.max()

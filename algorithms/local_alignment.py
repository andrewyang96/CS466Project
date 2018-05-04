"""Global alignment using the Needleman–Wunsch algorithm.

Source: http://schulte.faculty.unlv.edu//BIO480/SequenceAlign.pdf
"""

import numpy as np

NONE, LEFT, UP, DIAG = range(4)

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap_start_penalty=-2, gap_extend_penalty=-1):
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
    max_i, max_j, max_score = 0, 0, 0

    pointer[0, 0] = NONE
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

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

            if scores[i, j] <= 0:
                scores[i, j] = 0
                pointer[i, j] = NONE
                gap_i[i] = False
                gap_j[j] = False
            elif scores[i, j] > max_score:
                max_i = i
                max_j = j
                max_score = scores[i, j]

    aligned1 = []
    aligned2 = []
    i, j = max_i, max_j
    while scores[i, j] > 0:
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
    #print(scores)
    #print(pointer)
    #print(max_i, max_j, max_score)
    return ''.join(aligned1[::-1]), ''.join(aligned2[::-1]), scores.max()

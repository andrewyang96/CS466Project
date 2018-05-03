from algorithms.global_alignment import needleman_wunsch
from algorithms.local_alignment import smith_waterman
from algorithms.semiglobal_alignment import glocal_alignment

import pickle
import sys
import subprocess
import time

sys.setrecursionlimit(25000)

if __name__ == '__main__':
    with open('mutatedSeq.txt', 'rb') as f:
        seqs = pickle.load(f)

    print('Testing 20 global')
    global_results = []
    start = time.time()
    for i in range(20):
        global_results.append(needleman_wunsch(seqs[i][0], seqs[i][1]))
    end = time.time()
    global_time_elapsed = end - start

    print('Testing 20 local')
    local_results = []
    start = time.time()
    for i in range(20):
        local_results.append(smith_waterman(seqs[i][0], seqs[i][1]))
    end = time.time()
    local_time_elapsed = end - start

    print('Testing 20 glocal')
    glocal_results = []
    start = time.time()
    for i in range(20):
        glocal_results.append(glocal_alignment(seqs[i][0], seqs[i][1]))
    end = time.time()
    glocal_time_elapsed = end - start

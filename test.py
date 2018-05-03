from algorithms.global_alignment import needleman_wunsch
from algorithms.local_alignment import smith_waterman

import pickle
import sys
import subprocess
import time

sys.setrecursionlimit(25000)

if __name__ == '__main__':
    h_commands = ['python', 'algorithms/Hirschberg.py', 'mutatedSeq.txt', '-2', '-1', '2', '>', 'HirschbergResutls.txt']
    with open('mutatedSeq.txt', 'rb') as f:
        seqs = pickle.load(f)

    print('Testing global')
    start = time.time()
    global_result = needleman_wunsch(seqs[0], seqs[1])
    end = time.time()
    global_time_elapsed = end - start

    print('Testing local')
    start = time.time()
    local_result = smith_waterman(seqs[0], seqs[1])
    end = time.time()
    local_time_elapsed = end - start

    print('Testing Hirschberg')
    start = time.time()
    subprocess.call(h_commands)
    end = time.time()
    h_time_elapsed = end - start

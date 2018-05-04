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
    with open('mutatedSeqIDs.txt', 'rb') as f2:
        ids = pickle.load(f2)

    print('Testing 20 global')
    #global_results = []
    #start = time.time()
    for i in range(20):
        start = time.time()
        print("Aligning sequences " + str(ids[i][0]) + ' and ' + str(ids[i][1]) + ':' )
        #global_results.append(needleman_wunsch(seqs[i][0], seqs[i][1]))
        s1, s2, score = needleman_wunsch(seqs[i][0], seqs[i][1])
        print(s1)
        print(s2)
        print('Optimal score is '+ str(score))
        end = time.time()
        global_time_elapsed = end - start
        print('Global time elapsed: ' + str(global_time_elapsed))
        print("#" * 40)

    print('\n')
    print('Testing 20 local')
    #local_results = []
    for i in range(20):
        start = time.time()
        print("Aligning sequences " + str(ids[i][0]) + ' and ' + str(ids[i][1]) + ':' )
        #local_results.append(smith_waterman(seqs[i][0], seqs[i][1]))
        s1, s2, score = smith_waterman(seqs[i][0], seqs[i][1])
        print(s1)
        print(s2)
        print('Optimal score is '+ str(score))
        end = time.time()
        local_time_elapsed = end - start
        print('Local time elapsed: ' + str(local_time_elapsed))
        print("#" * 40)

    print('\n')
    print('Testing 20 glocal')
    #glocal_results = []
    for i in range(20):
        start = time.time()
        print("Aligning sequences " + str(ids[i][0]) + ' and ' + str(ids[i][1]) + ':' )
        #glocal_results.append(glocal_alignment(seqs[i][0], seqs[i][1]))
        glocal_alignment(seqs[i][0], seqs[i][1])
        end = time.time()
        glocal_time_elapsed = end - start
        print('Glocal time elapsed: ' + str(glocal_time_elapsed))
        print("#" * 40)

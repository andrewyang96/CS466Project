# CS466Project
Final project for CS 466 SP18.

For the project report, see the CS466 Project Report PDF.

DNA sequence data from http://genome.crg.es/datasets/genomics96/


The main dataset is "mutatedSeq.txt" which is in the format of a list of tuples. Each tuple contains 2 mutated sequences from data we retrieved from the given (above) link.
Sequence IDs of each sequence found in the main dataset can be found in "mutatedSeqIDs.txt" in the same format.

The shortened dataset is "deletedSeq.txt" which is in the format of a list of tuples.
Each tuple contains 2 mutated sequences. This is a modified version of the data also retrieved from the given (above) link where one of the sequences in the tuples is much shorter (~50%) than its pair sequence.
We did this so that we can compare how each alignment works when 2 sequences vary greatly in length.
Sequence IDs of each sequence found in the shortened dataset can be found in "deletedSeqIDs.txt" in the same format.

You can run global, local, and semiglobal alignments by first modifying the parameters
of each function in the test.py and/or test2.py.

  Run in the terminal (for main dataset):
    python test.py  

  Run in the terminal (for shortened dataset):
    python test2.py

You can run Hirschberg algorithm by modifying the arguments you pass in when typing in the terminal.

To run Hirschberg algorithm:
  python algorithm/Hirschberg.py <dataset .txt file> <dataset IDs .txt file> <gap_penalty> <mismatch> <match>

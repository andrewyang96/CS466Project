## Semi-global alignment, aka glocal alignment, according to CS 466 slides ##

import numpy as np

## returns optimal score and the aligned sequences ##
def glocal_alignment(seq1, seq2, match=2, gap_pen=-1, mismatch_pen=-2):
	#initialize final sequences
	aligned1 = ""
	aligned2 = ""

	len1 = len(seq1) + 1
	len2 = len(seq2) + 1

	#make sure seq2 is smaller than seq1
	if len1 < len2: 
		seq1, seq2 = seq2, seq1
		len1, len2 = len2, len1

	#make score and traceback matrices, initialize to zeros
	scores = np.zeros((len2, len1))
	traceback = np.zeros((len2, len1))

	#fill out score and traceback matrices
	for i in range(1, len2):
		for j in range(1, len1):
			if seq1[j-1] == seq2[i-1]:
				diagScore = scores[i-1, j-1] + match
			else:
				diagScore = scores[i-1, j-1] + mismatch_pen

			leftScore = scores[i, j-1] + gap_pen
			upScore = scores[i-1, j] + gap_pen

			if diagScore > leftScore and diagScore > upScore:
				idx = 0
				scores[i,j] = diagScore
			elif leftScore > diagScore and leftScore > upScore:
				idx = 1
				scores[i,j] = leftScore
			else:
				idx = 2
				scores[i,j] = upScore

			traceback[i, j] = idx #0 is diag, 1 is left, 2 is up in the traceback matrix
	
	maxIdx = np.argmax(scores[-1]) #index of max score

	#figures out actual sequences
	while i > 0 and j > 0:
		if j > maxIdx:
			aligned1 = seq1[j-1] + aligned1
			aligned2 = "-" + aligned2
			j = j-1
		else:
			tbScore = traceback[i][j]
			if tbScore == 0: #go up diagonally
				aligned1 = seq1[j-1] + aligned1
				aligned2 = seq2[i-1] + aligned2
				i = i-1
				j = j-1
			elif tbScore == 1: #go left
				aligned1 = seq1[j-1] + aligned1
				aligned2 = "-" + aligned2
				j = j-1
			elif tbScore == 2: #go up
				aligned2 = seq2[i-1] + aligned2
				aligned1 = "-" + aligned1
				i = i-1

	#fills out the front of sequence after alignment is done
	while i > 0:
		aligned1 = "-" + aligned1
		aligned2 = seq2[i-1] + aligned2
		i = i-1
	while j > 0:
		aligned2 = "-" + aligned2
		aligned1 = seq1[j-1] + aligned1
		j = j-1

	print("Max score:", max(scores[-1]))
	print("Sequences: \n" + aligned1 + "\n" + aligned2+"\n")
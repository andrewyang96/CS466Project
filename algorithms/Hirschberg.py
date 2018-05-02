import sys
import pdb
import argparse
import numpy as np
import copy
import re
import unittest
from Bio import SeqIO
import time

"""
Source: http://en.wikipedia.org/wiki/Hirschberg's_algorithm, https://github.com/wuzhigang05/Dynamic-Programming-Linear-Space
"""

arguments = sys.argv


insertion    = -500
deletion     = -500
substitution = -4
match = 5

def parse_fasta(fasta_file):
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #print(len(seq_record))  #print sequence length
        return seq_record.seq

def lastLineAlign(x, y):
  """
  input:  two strings: x and y
  output: an array with a length of |y| that contains the score for the alignment
          between x and y
  """
  global insertion
  global deletion
  global substitution
  global match
  row = y
  column = x
  minLen = len(y)
  prev = [0 for i in range(minLen + 1)]
  current = [0 for i in range(minLen + 1)]


  for i in range(1, minLen + 1):
    prev[i] = prev[i-1] + insertion

  current[0] = 0
  for j in range(1, len(column) + 1):
    current[0] += deletion
    for i in range(1, minLen + 1):
      if row[i-1] == column[j-1]:
        try:
          current[i] = max(current[i-1] + insertion, prev[i-1] + match, prev[i] + deletion)
        except:
          pdb.set_trace()
      else:
        current[i] = max(current[i-1] + insertion, prev[i-1] + substitution, prev[i] + deletion)
    prev = copy.deepcopy(current)

  return current

def partitionY(scoreL, scoreR):
  max_index = 0
  max_sum = float('-Inf')
  for i, (l, r) in enumerate(zip(scoreL, scoreR[::-1])):
    # calculate the diagonal maximum index
    if sum([l, r]) > max_sum:
      max_sum = sum([l, r])
      max_index = i
  return max_index

def dynamicProgramming(x, y):
  global insertion
  global deletion
  global substitution
  global match
  # M records is the score array
  # Path stores the path information, inside of Path:
  # d denotes: diagnal
  # u denotes: up
  # l denotes: left
  M = np.zeros((len(x) + 1, len(y) + 1))
  Path = np.empty((len(x) + 1, len(y) + 1), dtype=object)

  for i in range(1, len(y) + 1):
    M[0][i] = M[0][i-1] + insertion
    Path[0][i] = "l"
  for j in range(1, len(x) + 1):
    M[j][0] = M[j-1][0] + deletion
    Path[j][0] = "u"

  for i in range(1, len(x) + 1):
    for j in range(1, len(y) + 1):
      if x[i-1] == y[j-1]:
        M[i][j] = max(M[i-1][j-1] + match, M[i-1][j] + insertion, M[i][j-1] + deletion)
        if M[i][j] == M[i-1][j-1] + match:
          Path[i][j] =  "d"
        elif M[i][j] == M[i-1][j] + insertion:
          Path[i][j] = "u"
        else:
          Path[i][j] = "l"
      else:
        M[i][j] = max(M[i-1][j-1] + substitution, M[i-1][j] + insertion, M[i][j-1] + deletion)
        if M[i][j] == M[i-1][j-1] + substitution:
          Path[i][j] =  "d"
        elif M[i][j] == M[i-1][j] + insertion:
          Path[i][j] = "u"
        else:
          Path[i][j] = "l"

  row = []
  column= []
  middle = []
  i = len(x)
  j = len(y)
  while Path[i][j]:
    if Path[i][j] == "d":
      row.insert(0, y[j-1])
      column.insert(0, x[i-1])
      if x[i-1] == y[j-1]:
        middle.insert(0, '|')
      else:
        middle.insert(0, ':')
      i -= 1
      j -= 1
    elif Path[i][j] == "u":
      row.insert(0, '-')
      column.insert(0, x[i-1])
      middle.insert(0, 'x')
      i -= 1
    elif Path[i][j] == "l":
      column.insert(0, '-')
      row.insert(0, y[j-1])
      middle.insert(0, 'x')
      j -= 1
  return row, column, middle


def Hirschberg(x, y):
  row = ""
  column = ""
  middle = ""
#  x is being row-wise iterated (out-most for loop)
#  y is being column-wise iterated (inner-most of the for loop)
  if len(x) == 0 or len(y) == 0:
    if len(x) == 0:
      column = '-' * len(y)
      row = y
      middle =  'x' * len(y)
    else:
      column += x
      row += '-' * len(x)
      middle =  'x' * len(x)
  elif len(x) == 1 or len(y) == 1:
    row, column, middle = dynamicProgramming(x, y)
    # concatenate into string
    row, column, middle = map(lambda x: "".join(x), [row, column, middle])
  else:
    xlen = len(x)
    xmid = xlen/2
    ylen = len(y)

    scoreL = lastLineAlign(x[:int(xmid)], y)
    scoreR = lastLineAlign(x[int(xmid):][::-1], y[::-1])
    ymid = partitionY(scoreL, scoreR)
    row_l, column_u, middle_l = Hirschberg(x[:int(xmid)], y[:int(ymid)])
    row_r, column_d, middle_r = Hirschberg(x[int(xmid):], y[int(ymid):])
    row = row_l + row_r
    column = column_u + column_d
    middle = middle_l + middle_r

  return row, column, middle

if __name__ == '__main__':
  start_time = time.time()
  seqstr1 = ''.join(parse_fasta(arguments[1]))
  seqstr2 = ''.join(parse_fasta(arguments[2]))
  difference = []
  for i, (x, y) in enumerate(zip([seqstr1], [seqstr2])):
    row, column, middle = Hirschberg(x, y)
    difference.append(middle)
    print(row)
    print(middle)
    print(column)
    print
  alignment = ''.join(difference)
  num_matches = alignment.count('|')
  num_mismatches = alignment.count(':')
  num_gaps = alignment.count('x')
  score = num_matches * match + num_gaps * -500 + num_mismatches * substitution
  print(score)
  print("--- %s seconds ---" % (time.time() - start_time))

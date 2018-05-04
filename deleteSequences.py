from Bio import SeqIO
import sys
import random
from datetime import datetime
import math
import pickle

arguments = sys.argv

def write_to_textfile(file_name, list):
    f=open(file_name, "wb")
    #f.write(str(list))
    pickle.dump(list, f)
    f.close()

def mutate_sequence(seq):
    #print(seq)
    length = len(seq)
    #print("original_length is " + str(length))
    num_mutations = int(round(length/10))
    #print(num_mutations)
    for i in range(num_mutations):
        #print("length is " + str(len(modified_seq)))
        mutationIdx = random.randint(0, (len(seq) - 1))
        if random.uniform(0,1) <= 0.5: # random probability of mutation
            while True:
                random_int = random.randint(0,3)
                if(seq[mutationIdx] != nucleotides[random_int]): #makes sure nucleotide mutates to a different one
                    #print("mutate: " + str(mutationIdx))
                    #print(seq[mutationIdx], nucleotides[random_int])
                    seq[mutationIdx] = nucleotides[random_int] # mutate
                    break
        else:
            #print("delete: " + str(mutationIdx))
            del seq[mutationIdx]
    return ''.join(seq);

def delete_sequence(seq):
    length = len(seq)
    num_deletions = int(round(3*length/4))

    for i in range(num_deletions):
        deletionIdx = random.randint(0, (len(seq)-1))
        if random.uniform(0,1) <= 0.7:
            del seq[deletionIdx]
        return ''.join(seq)



    
if __name__ == '__main__':
    nucleotides = ['A', 'T', 'C', 'G']
    fasta_file = arguments[1]
    seq_tuple_list  = []
    seq_tuple_id_list = []
    with open(fasta_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            #print(list(''.join(record.seq)))
            S2 = mutate_sequence(list(''.join(record.seq)))
            S3 = delete_sequence(list(''.join(record.seq)))
            seq_tuple_list.append((''.join(S2),''.join(S3)))
            seq_tuple_id_list.append((record.id + "m1", record.id + "m2"))
    write_to_textfile("deletedSeq.txt", seq_tuple_list)
    write_to_textfile("deletedSeqIDs.txt", seq_tuple_id_list)

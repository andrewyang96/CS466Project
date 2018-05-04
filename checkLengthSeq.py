import pickle

file = open('deletedSeq.txt', 'rb')
result = pickle.load(file)

file2 = open('mutatedSeq.txt', 'rb')
result2 = pickle.load(file2)

print("deletedSequences")
for items in result:
    print((len(items[0]), len(items[1])))

print("mutatedSequences")
for items in result2:
    print((len(items[0]), len(items[1])))

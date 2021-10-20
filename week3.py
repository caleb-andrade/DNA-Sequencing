# Week 3. Practical: implementing dynamic programming for edit distance
import time

def readGenome(filename):
    """
    Parses a DNA reference genome from a file in the FASTA format,
    Authors: Ben Langmead & Jacob Pritt.
    """
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
def readFastq(filename):
    """
    Parses read and quality strings from a FASTQ file with sequencing reads
    Authors: Ben Langmead & Jacob Pritt.
    """
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() #read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def editDistance(x, y):
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))

    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1

            D[i][j] = min (distHor, distVer, distDiag)

    return D[-1][-1]

alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8], [4, 0, 4, 2, 8],
         [2, 4, 0, 4, 8], [4, 2, 3, 0, 8],
         [8, 8, 8, 8, 8]]

def globalAlignment(x, y):
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))

    for i in range(1, len(x)+1):
        D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]
    for i in range(1, len(y)+1):
        D[0][i] = D[0][i-1] + score[-1][alphabet.index(y[i-1])]
        
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            D[i][j] = min (distHor, distVer, distDiag)

    return D[-1][-1]

def overlap(a, b, min_length=3):
    start = 0

    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1

def naive_overlap_map(reads, k):
    olaps = {}
    for a,b in permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps

# Testing edit distance
x = 'shake spea'
y = 'Shakespear'
print editDistance(x, y)

# Testing global Alignment
x = 'TACCAGATTCGA'
y = 'TACCAGATTCGC'
print globalAlignment(x, y)

#Testing overlap
print overlap('TTACGT', 'CGTACCG')
from itertools import permutations
print list(permutations([1,2,3], 2))

reads = ['ACGGATTATC', 'GATCAAGT', 'TTCACGGA']
print naive_overlap_map(reads, 3)


def editDistanceAppMat(p, t):
    D = []
    # initialize matrix    
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    # fill first column of matrix
    for i in range(len(p)+1):
        D[i][0] = i
    # fill first row of 0's
    for i in range(len(t)+1):
        D[0][i] = 0

    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1

            D[i][j] = min (distHor, distVer, distDiag)
    
    # search for the minimum value in last row
    best = len(p)
    for col in range(len(t) + 1):
        if D[len(p)][col] < best:
            best = D[len(p)][col]
            
    return best
    
def overlapMap(sequences, k):
    
    kmers_dict = {}
    tic = time.clock()
    # construct the dictionary
    for read in sequences:
        for i in range(1+len(read)-k):
            kmer = read[i:k+i]
            if kmers_dict.has_key(kmer):
                kmers_dict[kmer].add(read)
            else:
                kmers_dict[kmer] = set([read])
    
    # check matches with suffix of every read
    overlaps = 0
    hits = 0
    for read in sequences:
        hit = False
        for test_read in kmers_dict[read[-k:]]:
            if test_read != read:
                if overlap(read, test_read, min_length = k) > 0:
                    overlaps += 1
                    hit = True
        if hit:
            hits += 1
            
    toc = time.clock()
            
    return overlaps, hits, toc-tic 
    
text = readGenome('chr1.GRCh38.excerpt.fasta')
sequences = readFastq('ERR266411_1.for_asm.fastq')[0]

print "\nQuestion 1"
pattern1 = 'GCTGATCGATCGTACG'
print "Edit distance: ", editDistanceAppMat(pattern1, text)

print "\nQuestion 2"
pattern2 = 'GATTTACCAGATTGAG'
print "Edit distance: ", editDistanceAppMat(pattern2, text)

#sequences = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
#sequences = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
answer = overlapMap(sequences, 30)
print "\nQuestion 3"
print "Number of overlaps: ", answer[0]
print "Running time: ", answer[2]

print "\nQuestion 4"
print "Number of hits: ", answer[1]



                                           

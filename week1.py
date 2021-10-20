# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 15:03:39 2015

@author: Caleb Andrade
"""

import matplotlib.pyplot as plt
import collections


def naive(p, t):
    """
    Naive algorithm for exact matching,
    Authors: Ben Langmead & Jacob Pritt.
    """
    ocurrences = []
    # loop over all alignments
    for i in range(len(t) - len(p) + 1):
        match = True
        # loop over characters
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            ocurrences.append(i)
    return ocurrences
    
def reverseComplement(s):
    """
    Takes DNA string and returns its reverse complement.
    Authors: Ben Langmead & Jacob Pritt.
    """
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
    
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

def naive_with_rc(p, t):
    """
    Naive algorithm for exact matching, considering reverse complement.
    """
    ocurrences = []
    def compare(p, t, i):
        """
        Nested function to loop and compare a patterns characters.
        """
        for j in range(len(p)):
            if t[i+j] != p[j]:
                return False
        return True
        
    # loop over all alignments
    for i in range(len(t) - len(p) + 1):
        if compare(p, t, i):
            ocurrences.append(i)
            continue
        rc = reverseComplement(p)
        if compare(rc, t, i):
            ocurrences.append(i)
    return ocurrences
    
def naive_2mm(p, t):
    """
    Naive algorithm for matching allowing 2 mismatches
    Authors: Ben Langmead & Jacob Pritt.
    """
    ocurrences = []
    # loop over all alignments
    for i in range(len(t) - len(p) + 1):
        match = True
        # loop over characters
        flag = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                flag += 1
                if flag > 2:
                    match = False
                    break
        if match:
            ocurrences.append(i)
    return ocurrences
    
def findQualityByPos(quals):
    """
    """
    qp = [0]*100
    totals = [0]*100
    for qual in quals:
        for i in range(len(qual)):
            qp[i] += phred33ToQ(qual[i])
            totals[i] += 1
                
    for i in range(len(qp)):
        qp[i] /= float(totals[i])
            
    return qp

def phred33ToQ(qual):
    return ord(qual) - 33
    
def createHist(qualities):
    hist = [0]*50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

# testing naive
#p = 'ACG'
#t = 'ACGACTGGTCACGTACTGACGT'
#print naive(p, t)
 
# testing reverse Complement
#rc = reverseComplement(p)
#print rc
#print naive(rc, t)

# testing naive with rc
#print naive_with_rc(p, t)

# testing naive_2mm
#p = 'ACTTTA'
#t = 'ACTTACTTGATAAAGT'
#print naive_2mm(p, t)


# testing readGenome
lvg = readGenome('lambda_virus.fa')

# Question 1
print "\nQuestion 1"
print len(naive_with_rc('AGGT', lvg))

# Question 2
print "\nQuestion 2"
print len(naive_with_rc('TTAA', lvg))

# Question 3
print "\nQuestion 3"
print naive_with_rc('ACTAAGT', lvg)[0]

# Question 4
print "\nQuestion 4"
print naive_with_rc('AGTCGA', lvg)[0]

# Question 5
print "\nQuestion 5"
print len(naive_2mm('TTCAAGCC', lvg)) 

# Question 6
print "\nQuestion 6"
print naive_2mm('AGGAGGTT', lvg)[0]

# Question 7
print "\nQuestion 7"
seqs, quals = readFastq('ERR037900_1.first1000.fastq')

h = createHist(quals)
plt.bar(range(len(h)), h)
plt.show()

qp = findQualityByPos(quals)
for i in range (len(qp)):
    print i, " ", qp[i]
plt.plot(range(len(qp)), qp)
plt.show()

count = collections.Counter()
for seq in seqs:
    count.update(seq)
print(count)



                
    
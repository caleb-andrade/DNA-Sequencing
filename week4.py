# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:30:08 2015

@author: Ben Langmead and Jacob Pritt, unless stated otherwise.
"""

import itertools
import time

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

def overlap(a, b, min_length=3):
    """
    Return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0. 
    """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def scs(ss):
    """ 
    Returns shortest common superstring of given
    strings, which must be the same length. 
    """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            optimal_set = [sup]
        elif len(sup) == len(shortest_sup):
            optimal_set.append(sup)
    
    return shortest_sup, optimal_set  # return shortest
    
def scs_mod(reads, k):
    """ 
    Returns shortest common superstring of given
    strings, which must be the same length.
    """
    # Create a dictionary. Every read is a key and its value is a list
    # of tuples. Each tuple represents an overlapping read and the size
    # of the overlapping (size, read).
    tic = time.clock()    
    olen_dict = {}
    indices = set(i for i in range(len(reads)))
    for i in range(len(reads)):
        olen_dict[i] = []
        for j in range(len(reads)):
            if reads[i] != reads[j]:
                olen = overlap(reads[i], reads[j], min_length = k)
                if olen != 0:
                    olen_dict[i].append((olen, j))
                    indices.discard(j)
     
    # now, create an ordered list (max overlapping) to construct superstring.    
    idx = list(indices)[0]
    ordlist = [reads[idx]]
    for i in range(len(olen_dict)):
        temp = olen_dict[idx]
        temp.sort(reverse = True)  
        if len(temp) == 0:
            break
        size_read = temp[0]
        ordlist.append(reads[size_read[1]])
        idx = size_read[1]
 
    sup = ordlist[0]  # superstring starts as first string
    for i in range(len(ordlist)-1):
        # overlap adjacent strings A and B in the permutation
        olen = overlap(ordlist[i], ordlist[i+1])
        # add non-overlapping portion of B to superstring
        sup += ordlist[i+1][olen:]
    
    toc = time.clock()
    print "Running time: ", toc - tic
        
    return sup  # return shortest
    
print "\nQuestion 1: "
ss = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print len(scs(ss)[0])

print "\nQuestion 2: "
print scs(ss)[1]

print "\nQuestion 3: "
vr = readFastq('ads1_week4_reads.fq')[0]
#print "Number of reads: ", len(vr)
scs = scs_mod(vr, 30)
print "Length of scs: ", len(scs)
countA = 0
countT = 0
for char in scs:
    if char == 'A':
        countA += 1
    if char == 'T':
        countT += 1
print "Number of A's in scs: ", countA

print "\nQuestion 4: "
print "Number of T's in scs: ", countT

print "\nQuestion 5: Measles"
    
    
    
    

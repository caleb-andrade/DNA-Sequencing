# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:51:08 2015

@author: Benjamin Langmead & Jacob Pritt, unless stated otherwise.
"""
import bm_preproc as bmp

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

def boyer_moore(p, p_bm, t):
    """
    Do Boyer-Moore matching.
    """
    num_char = 0 # count for number of character comparisons
    num_alig = 0 # count for number of alignments tried
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        num_alig += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            num_char += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    
    print "\nBoyer-Moore algorithm"
    print "Number of character comparisons: ", num_char
    print "Number of alignments tried: ", num_alig

    return occurrences
    

# INDEXING

import bisect

class Index(object):
    def __init__(self, t, k):
        """
        Create index from all substrings of size 'length'
        """
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()
    
    def query(self, p):
        """
        Return index hits for first k-mer of P
        """
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        
def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k:i+len(p)]:
            offsets.append(i)
    return offsets
    
# APPROXIMATE MATCHING

def approximate_match(p, t, n):
    """
    Approximate matching using boyer-moore. Pigeon-hole principle.
    """
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = bmp.BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)
    
def naive(p, t):
    """
    Naive algorithm for exact matching,
    Authors: Ben Langmead & Jacob Pritt.
    """
    num_char = 0 # count for number of character comparisons
    num_alig = 0 # count for number of alignments tried
    ocurrences = []
    # loop over all alignments
    for i in range(len(t) - len(p) + 1):
        num_alig += 1
        match = True
        # loop over characters
        for j in range(len(p)):
            num_char += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            ocurrences.append(i)
            
    print "\nNaive algorithm"
    print "Number of character comparisons: ", num_char
    print "Number of alignments tried: ", num_alig
    
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

def app_match_idx(p, t, n):
    """
    Approximate matching using index. Pigeon-hole principle.
    'n' stands for the number of allowed mismatches.
    """
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    index = Index(t, 8) # 8-mer index
    hits = 0 # keep track of hits before verification
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = queryIndex(p[start:end], t, index)
        hits += len(matches)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    if mismatches > n:
                        break
                    mismatches += 1
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
                
    print "Hits: ", hits
    
    return list(all_matches)

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        
def app_match_subseq_idx(p, t, n):
    """
    Approximate matching using sub-seq-index. 
    'n' stands for the number of allowed mismatches.
    """
    all_matches = set()
    ind = SubseqIndex(t, 8, 3)
    hits = 0
    for i in range(3):
        matches = ind.query(p[i:])
        hits += len(matches)
        # verify if there are matchings with p with at most 2 mismatches
        for m in matches:            
            if not m + len(p) - i < len(t):
                continue
            mismatches = 0
            
            for j in range(len(p)):
                if not p[j] == t[m + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - i)
    print "Hits: ", hits
    return list(all_matches)

#**********************************************************************  
print "\nQuestion 1, 2, 3"
t = readGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p_bm = bmp.BoyerMoore(p, alphabet='ACGT')
print len(naive(p, t))
print boyer_moore(p, p_bm, t)
#**********************************************************************  
print "\nQuestion 4, 5"
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
matchings = app_match_idx(p, t, 2)
#matchings = naive_2mm(p, t)
matchings.sort()
print "Matchings: ", matchings
print "Number of matchings: ", len(matchings)
#**********************************************************************  
print "\nQuestion 6"
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
matchings = app_match_subseq_idx(p, t, 2)
#print len(naive(p, t))
print "Matchings: ", matchings
print "Number of matchings: ", len(matchings)
#**********************************************************************  
"""
# Examples of bad_character and good_suffix methods
p = 'TCAA'
p_bm = bmp.BoyerMoore(p, alphabet = 'ACGT')
print p_bm.bad_character_rule(2, 'T')

p = 'ACTA'
p_bm = bmp.BoyerMoore(p, alphabet = 'ACGT')
print p_bm.good_suffix_rule(0)

p = 'ACAC'
p_bm = bmp.BoyerMoore(p, alphabet='ACGT')
print p_bm.match_skip()

# Example of Boyer-Moore method
t = 'GCTAGCTCTACGAGTCTA'
p = 'TCTA'
p_bm = bmp.BoyerMoore(p, alphabet='ACGT')
print boyer_moore(p, p_bm, t)

# Example of index querying
t = 'ACTTGGAGATCTTTGAGGCTAGGTATTCGGGATCGAAGCTCATTTCGGGGATCGATTACGATATGGTGGGTATTCGGGA'
p = 'GGTATTCGGGA'
index = Index(t, 4)
print "\nIndex querying example"
print(queryIndex(p, t, index))

# Example of approximate matching
p = 'AACTTG'
t = 'CACTTAATTTG'
print "\nApproximate matching example"
print(approximate_match(p, t, 2))

# Examples of subsequence indexing
ind = SubseqIndex('GGGGATATATTCTGTCG', 3, 2)
print (ind.index)

p = 'TTATATCTC'
print(ind.query(p[0:]))
print(ind.query(p[1:]))
print(ind.query(p[2:]))
print(ind.query(p[4:]))
print(ind.query(p[5:]))

# Example of subsequence approximate matching
t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
print app_match_subseq_idx(p, t, 2)
"""


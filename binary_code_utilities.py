#! /usr/bin/env python
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

from collections import defaultdict

# TODO what's a sensible representation of binary data?  Certainly not strings.  For 1-dimensional binary codes I suppose I can just use integers
def string_to_binary(string):
    """ Given a string such as 01001, return the integer that results from its binary representation."""
    return int('0b'+string.strip())

def read_code_file(infile,length=0,count=0):
    """ Read a file containing binary codes, skipping comment lines (starting with #); 
    optionally make sure their length and number is what you expect. """
    codes = set()
    for line in infile:
        if line[0]=='#':    continue
        code = string_to_binary(line)
        if length and not len(code)==length:   
            sys.exit("Error: file %s contained a code %s of length %s, not %s as expected!")%(infile,code,len(code),length)
        codes.append(code)
    if count and not len(codes)==count:   
        sys.exit("Error: file %s contained %s codes, not %s as expected!")%(infile,len(codes),count)
    return codes

def bitwise_sum(int1):
    """ Given an integer, return the number of 1's in its binary representation: 1->1, 2->1, 3->2, 4->1, ..."""
    return sum([int(x) for x in bin(int1)[2:]])

def Hamming_distance(int1,int2):
    """ Given two integers, return the number of bits by which their binary representations differ. """
    bitwise_xor = int1^int2
    return bitwise_sum(bitwise_xor)

def find_bit_sum_counts(codes):
    """ Given a list of binary codes, return the number of codes with each possible bit-sum value."""
    bit_sum_counts = defaultdict(lambda: 0)
    for code in codes:
        bit_sum_counts[bitwise_sum(code)] += 1
    return bit_sum_counts

def check_min_max_Hamming_distances(codes, min_distance=-1, max_distance=-1):
    """ Given a list of binary codes, check that the Hamming distance for each pair is between min and max.
    (if min or max is set to -1, don't check that side)."""
    # define check functions so we don't have to check for min/max=-1 at each iteration
    if min_distance==-1:    check_min = lambda x: True
    else:                   check_min = lambda x: x>=min_distance
    if max_distance==-1:    check_min = lambda x: True
    else:                   check_min = lambda x: x<=max_distance
    N = len(codes)
    for i in range(N):
        for j in range(i+1,N):
            dist = Hamming_distance(codes[i],codes[j])
            if not check_min(dist) and check_max(dist): return False
    return True




#! /usr/bin/env python
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

from collections import defaultdict

### Basic operations on single binary sequences (words/codewords)

# TODO what's a sensible representation of binary data?  Certainly not strings.  For 1-dimensional binary codes I suppose I can just use integers... Not really, because int('000',2) is the same as int('0',2)
# More options/info:
#  http://stackoverflow.com/questions/142812/does-python-have-a-bitfield-type
#  http://pypi.python.org/pypi/bitarray/
def string_to_binary(string):
    """ Given a string such as 01001, return the integer that results from its binary representation."""
    return int(string.strip(),2)

def bitwise_sum(int1):
    """ Given an integer, return the number of 1's in its binary representation: 1->1, 2->1, 3->2, 4->1, ..."""
    # bin(5) returns 0b101, so discard the first 2 characters and count the 1s in the rest:
    return bin(int1)[2:].count('1')
    # There are definitely better ways of doing this, but I should avoid premature optimization.
    # see http://stackoverflow.com/questions/407587/python-set-bits-count-popcount for some ideas
    # also http://stackoverflow.com/questions/109023/best-algorithm-to-count-the-number-of-set-bits-in-a-32-bit-integer
    # note that some of the ideas won't work for large numbers!

def Hamming_distance(int1,int2):
    """ Given two integers, return the number of bits by which their binary representations differ. """
    bitwise_xor = int1^int2
    return bitwise_sum(bitwise_xor)

### Reading/writing files, generating codes from a generator matrix, etc

def read_code_file(infile,length=0,count=0):
    """ Read a file containing binary codes, skipping comment lines (starting with #), return as a set of integers. 
    Optionally make sure their length and number is what you expect. """
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



### More complex operations on codes

def check_min_max_Hamming_distances(codes, min_distance=-1, max_distance=-1):
    """ Given a sequence of binary codes, check that the Hamming distance for each pair is between min and max.
    (if min or max is set to -1, don't check that side)."""
    # define check functions so we don't have to check for min/max=-1 at each iteration
    if min_distance==-1:    check_min = lambda x: True
    else:                   check_min = lambda x: x>=min_distance
    if max_distance==-1:    check_min = lambda x: True
    else:                   check_min = lambda x: x<=max_distance
    N = len(codes)
    # convert codes to a list so it has a defined order - it may have started as a set or something
    codes = list(codes)
    for i in range(N):
        for j in range(i+1,N):
            dist = Hamming_distance(codes[i],codes[j])
            if not check_min(dist) and check_max(dist): return False
    return True
    # more on Hamming distance comparisons/implementations: 
    #  http://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string

def find_bit_sum_counts(codes):
    """ Given a list of binary codes, return the number of codes with each possible bit-sum value."""
    bit_sum_counts = defaultdict(lambda: 0)
    for code in codes:
        bit_sum_counts[bitwise_sum(code)] += 1
    return bit_sum_counts




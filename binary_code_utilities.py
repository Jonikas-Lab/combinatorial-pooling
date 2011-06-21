#! /usr/bin/env python
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

from collections import defaultdict
import bitstring

class BinaryCodeError(Exception):
    """ Exceptions in the binary_code_utilities module."""
    pass

### Binary string representations and basic functions
# what's a sensible representation of binary data?  Certainly not strings.  For 1-dimensional binary codes I suppose I can just use integers... Not really, because int('000',2) is the same as int('0',2)
#  http://stackoverflow.com/questions/142812/does-python-have-a-bitfield-type
# should support basic bitwise operators: &, |, ^, ~

class Binary_codeword:
    """ A binary codeword like '001' or '0110101010101'. Nothing complicated. """
    # implemented with bistring: http://pypi.python.org/pypi/bitstring/2.2.0  http://code.google.com/p/python-bitstring/
    # more convenient/logical than bitarray; has more documentation/support; may be slower but this shouldn't matter much.

    def __init__(self,val,length=0,check_length=False):
        """ Generate the self.codeword binary word based on val; pad with 0s on the left to desired length if specified.
        If check_length is True, instead of padding make sure the length is as specified, raise BinaryCodeError if not.
        If val is a 0/1 string, strip spaces/newlines and convert straight to a bit-string.
        If val is a list of 0/1 or True/False values, 
        If val is an int, the bitstring should be as if the builtin bin function was used and the initial 0b stripped.
        How other argument types will behave is not guaranteed - depends on the package used for the representation."""
        # for binary strings just specify it's bin: '110' and '0b110' and '  110\n' and '11 0' all give the same result
        if isinstance(val,str):   self.codeword = bitstring.BitArray(bin=val)
        # for ints, use the uint method if we know the length, otherwise have to convert to string
        elif isinstance(val,int):     
            if length:  self.codeword = bitstring.BitArray(uint=val,length=length)
            else:       self.codeword = bitstring.BitArray(bin(val))
        # lists of 0/1 or True/False values natively work as expected; I don't know or care what happens with other types.
        else:                       self.codeword = bitstring.BitArray(val)   
        self.codeword = bitstring.BitArray(val)
        # pad to given length or check the length if necessary
        if length and not length==len(self):
            if not check_length:
                self.pad(length,0)
            else:
                if not self.check_length(length):
                    raise BinaryCodeError("The created binary codeword didn't match the expected length!")

    def pad(self,length,value=0):
        """ Pad on the left to the given length (with 0 by default). """
        length_diff = length - len(self)
        # if the length is already correct, do nothing; if it's too high, complain
        if length_diff==0:      pass
        elif length_diff < 0:   raise BinaryCodeError("Can't pad the codeword to a length lower than its current length!")
        else:                   self.codeword.prepend('0b' + length_diff * str(int(value)))

    def check_length(self,length):
        """ Pad on the left to the given length (with 0 by default). """
        return length==len(self.codeword)

    def weight(self):
        """ Return the number of 1's in the codeword."""
        return self.codeword.count(1)

    def string(self):
        """ Return a plain 0/1 string representation. """
        return self.codeword.bin[2:]

    def __len__(self):
        """ Return the length of the codeword."""
        return self.codeword.length()

def Hamming_distance(val1,val2):
    """ Given two binary strings, return the number of bits by which their binary representations differ. """
    bitwise_xor = val1^val2
    return bitwise_xor.weight()


### Binary code (set of binary strings) representation

class Binary_code:
    # TODO implement the rest of this!  Just started

    def __init__(self,length):
        self.length = length
        self.codewords = set()

    def add(self,val):
        new_word = binary_codeword(val)
        if not len(new_word)==self.length: 
            raise BinaryCodeError("Trying to add a codeword of the wrong length to a binary code!")
        self.codewords.add(new_word)

    def size(self):
        return len(self.codewords)

    # TODO some of the more complex functions from below should actually be in this class


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




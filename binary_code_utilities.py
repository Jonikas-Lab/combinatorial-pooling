#! /usr/bin/env python2
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

from collections import defaultdict
from itertools import combinations
from numpy import array, dot
import bitstring

class BinaryCodeError(Exception):
    """ Exceptions in the binary_code_utilities module."""
    pass

### Binary string representations and basic functions
# what's a sensible representation of binary data?  Certainly not strings.  For 1-dimensional binary strings I suppose I can just use integers... Not really, because int('000',2) is the same as int('0',2)
#  http://stackoverflow.com/questions/142812/does-python-have-a-bitfield-type
# should support basic bitwise operators: &, |, ^, ~

class Binary_codeword:
    """ A binary string like '01101'. '001' and '01' are distinct. Supports |, &, ^, ~ bitwise operators.
    Not just a binary representation of an integer: '001' and '01' are distinct. """
    # implemented with bistring: http://pypi.python.org/pypi/bitstring/2.2.0  http://code.google.com/p/python-bitstring/
    # more convenient/logical than bitarray; has more documentation/support; may be slower but this shouldn't matter much.

    def __init__(self,val,length=0,check_length=False):
        """ Generate the self.codeword binary word based on val; pad with 0s on the left to desired length if specified.
        If check_length is True, instead of padding make sure the length is as specified, raise BinaryCodeError if not.
        If val is a 0/1 string, strip spaces/newlines and convert straight to a binary string.
        If val is an int, act as if the builtin bin() function was used to convert to a string (stripping initial 0b).
        If val is a list of 0/1 or True/False values, treat it as the corresponding 0/1 string.
        If val is a Binary_codeword instance, just copy it.
        How other argument types will behave is not guaranteed - depends on the package used for the representation."""
        # for binary strings just specify it's bin: '110' and '0b110' and '  110\n' and '11 0' all give the same result
        if isinstance(val,str):   
            self.codeword = bitstring.BitArray(bin=val)
        # for ints, use the uint method if we know the length, otherwise have to convert to string
        elif isinstance(val,int):
            if length:  self.codeword = bitstring.BitArray(uint=val,length=length)
            else:       self.codeword = bitstring.BitArray(bin(val))
        # for Binary_codeword instances, make a new one with the same bitstring
        elif isinstance(val,Binary_codeword):
            self.codeword = bitstring.BitArray(val.codeword)
        # lists of 0/1 or True/False values natively work as expected; I don't know or care what happens with other types.
        else:                       self.codeword = bitstring.BitArray(val)   
        # pad to given length or check the length if necessary
        if length and not length==len(self):
            # if only checking the length, raise an error
            if check_length:
                raise BinaryCodeError("The created binary codeword didn't match the expected length!")
            # otherwise pad with 0s on the left to the given length (does nothing if length < len(self)
            self.codeword.prepend((length-len(self))*[0])
            # unless it's already too long, then raise an error
            if length-len(self) < 0:
                raise BinaryCodeError("Can't pad the codeword %s to length %s, "%(self.codeword,length)
                                      + "since that's lower than its current length!")

    def weight(self):
        """ Return the number of 1's in the codeword (i.e. the Hamming weight or bitwise sum)."""
        return self.codeword.count(1)

    def string(self):
        """ Return a plain 0/1 string representation. """
        return self.codeword.bin[2:]

    def list(self):
        """ Return a representation as a list of ints (0 or 1). """
        return [int(x) for x in list(self.codeword)]

    def __len__(self):
        """ Return the length of the codeword."""
        return self.codeword.length

    # Bitwise and, or, xor, not operators - just passed on to the underlying bitstring objects.
    def __and__(self,other):    return Binary_codeword(self.codeword & other.codeword)
    def __or__(self,other):     return Binary_codeword(self.codeword | other.codeword)
    def __xor__(self,other):    return Binary_codeword(self.codeword ^ other.codeword)
    def __invert__(self):       return Binary_codeword(~self.codeword)

    def __cmp__(self,other):
        """ Comparison/sorting based on string representation: Two instances with the same bitstring should be equal."""
        # NOTE: comparison and hashing are related and need to match!
        #   By default, they're based on identity: hash(x) is id(x), and x==y iff id(x)==id(y), i.e. if 
        #     x and y are the SAME object, not just different objects with identical contents. 
        #   So unless I overload the operators, a set can contain two distinct Binary_code('111') instances! BAD.
        #   If I implement __cmp__ but not __hash__, the objects are considered unhashable, because otherwise
        #     there can be cases where x==y but hash(x)!=hash(y), which is BAD for hashing.  
        #   See http://docs.python.org/reference/datamodel.html#object.__hash__ for more on this.
        #   TODO really, in order to make this absolutely right, I should make sure the bitstrings are immutable...
        # TODO is string comparison really what I want here?  How about when the lengths are different? Should bitstrings
        #   with different lengths even be comparable?  I suppose they should just so I can sort stuff and get 
        #   a consistent result. Possibly just sorting by length first would be better, but it doesn't matter much.
        #   As long as identity works correctly and there's SOME sensible sorting, we're fine.
        return cmp(self.string(),other.string())

    def __hash__(self):
        """ Hashing based on the string representation."""
        # NOTE: comparison and hashing are related and need to match!  See notes in __cmp__
        return hash(self.string())

    # This is exactly as it should be: __repr__ gives the full evaluatable definition, __str__ gives something readable.
    def __repr__(self):     return "Binary_codeword('%s')"%self.string()
    def __str__(self):      return self.string()

def Hamming_distance(val1,val2):
    """ Given two binary strings, return the number of bits by which their binary representations differ. """
    bitwise_xor = val1^val2
    return bitwise_xor.weight()


### Binary code (set of binary strings) representation

# TODO figure out a naming that won't confuse people!  (Or me!)  Mathematically a "code" is a set of codewords (or a method of encoding things), but IRL a "code" can be either a method of encoding things or just a string, so code/codeword is confusing and even I'm using them wrong!

class Binary_code:
    """ Essentially a set of Binary_codeword objects, all of the same length."""

    def __init__(self,length,val=[],method='list',expected_count=0):
        """ Initialize with the given length, and add all elements of values to the set of codewords (default empty)."""
        try: 
            self.length = int(length)
        except ValueError:  
            raise BinaryCodeError("The length argument to Binary_code must be an int or possible to cast to an int!")
        self.method = method
        self.codewords = set()
        if method=='list':
            for x in val: self.add(x)
        elif method=='listfile':    
            self.read_code_from_file(val)
        elif method=='matrix':    
            self.get_code_from_generator_matrix(generator_matrix=val)
        elif method=='matrixfile':    
            self.get_code_from_generator_matrix(generator_file=val)
        else:
            raise BinaryCodeError("method %s passed to the Binary_code initializer not recognized! "%method 
                                  + "(allowed methods are list, listfile, matrix, matrixfile)")
        if expected_count and not self.size()==expected_count:
            raise BinaryCodeError("Binary_code initializer gave %s codewords, "%self.size() 
                                  + "not %s as expected!"%expected_count)

    def add(self,val):
        """ Add Binary_code(val) codeword to the code, checking for correct length."""
        # it's all right if val is a Binary_codeword already, that works too - similar to sets
        self.codewords.add(Binary_codeword(val,length=self.length,check_length=True))

    def remove(self,val):
        """ Remove Binary_code(val) codeword from the code; fail if val wasn't in the code, or is the wrong length."""
        try:
            self.codewords.remove(Binary_codeword(val,length=self.length,check_length=True))
        # trying to remove an element from a set where it wasn't present raises KeyError - we want a similar behavior.
        except KeyError:
            raise BinaryCodeError("Codeword %s cannot be removed from code because it wasn't present!"%val)

    def remove_all_zero_codeword(self):
        """ If the code contains the all-zero codeword, remove it and return 1; otherwise return 0; never fail. """
        try:
            self.codewords.remove(Binary_codeword('0'*self.length,length=self.length))
            return 1
        except KeyError:
            return 0

    def size(self):
        """ Return the number of codewords currently in the code."""
        return len(self.codewords)

    def read_code_from_file(self,infile,expected_count=0):
        """ Populate the code with codewords read from a plaintext file of 0/1 strings (one per line).
        Skip comment lines (starting with #).  Optionally make sure the codeword count is as expected. """
        read_count = 0
        for line in open(infile):
            if not line.startswith('#'):
                self.add(line)
                read_count += 1
        if expected_count and not read_count==expected_count:   
            raise BinaryCodeError("File %s contained %s codewords, "%(infile,read_count)
                                  + "not %s as expected!"%expected_count)

    def write_code_to_file(self,outfile):
        """ Write all the codewords (in arbitrary order) to a plaintext file of 0/1 strings (one per line). """
        for codeword in self.codewords:
            outfile.write(codeword.string()+'\n')
        # TODO actually write a code to a file (preferably one we'll use), test writing and reading!

    def get_code_from_generator_matrix(self,generator_file=None,generator_matrix=None,expected_count=0):
        """ Given either a generator matrix (as a numpy array) or a plaintext file containing the matrix, 
        add all codewords generated by that matrix to the current code.  Check if the count is as expected, if given.
        More information on how a generator matrix works: http://en.wikipedia.org/wiki/Generator_matrix."""
        # make sure we didn't get a file AND an arra; read the matrix from the file if it was provided as a file
        if generator_file and generator_matrix:
            raise BinaryCodeError("Provide either a generator_file or a generator_matrix, not both!")
        if generator_file:
            generator_matrix = array([[int(x) for x in line.strip()] for line in open(generator_file)])
        # the encoded word length and the codeword length (i.e. the n and k from the [n,k,d] description)
        #  can be inferred from the generator matrix size; check if they match self.length and expected_count (if given)
        input_code_length, output_code_length = generator_matrix.shape
        codeword_count = 2**input_code_length
        if not self.length==output_code_length:
            raise BinaryCodeError("Trying to use a generator matrix of shape (%s,%s) "%generator_matrix.shape
                                  + "to generate codewords of length %s - the sizes don't match up!"%self.length)
        if expected_count and not expected_count==codeword_count:
            raise BinaryCodeError("The number of codewords generated by the given matrix will be %s, "%codeword_count
                                  + "not %s as expected!"%expected_count)
        # generate all possible input codewords of length k, and convert them to codewords using matrix dot-multiplication
        #   (need to do a %2 at the end because binary addition in coding theory is basically xor - I think...)
        for x in range(codeword_count):
            x = Binary_codeword(x,length=input_code_length).list()
            self.add(dot(x,generator_matrix)%2)

    # TODO I could make one or both of these be the initialization signature instead, but who cares
    # TODO could add the minimum Hamming distance to this?
    def __str__(self):  return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())
    def __repr__(self): return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())

    def find_Hamming_distance_range(self):
        """ Return a tuple containing the lowest and highest Hamming distance between all codeword pairs."""
        # set start lov/high values to impossible extremes to make sure they get overwritten
        low = self.length+1
        high = 0
        for x,y in combinations(self.codewords,2):
            dist = Hamming_distance(x,y)
            low = min(low,dist)
            high = max(high,dist)
        return low,high
        # more on Hamming distance comparisons/implementations: http://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string

    def find_bit_sum_counts(self):
        """ Return the number of codewords with each possible bit-sum value (weight).
        The return value is a sorted list of tuples, for readability, but convertion into a dictionary is trivial."""
        bit_sum_counts = defaultdict(lambda: 0)
        for codeword in self.codewords:
            bit_sum_counts[codeword.weight()] += 1
        return sorted(bit_sum_counts.items())

    def total_bit_sum(self):
        """ Return the total sum of bits in all the codewords."""
        return sum([codeword.weight() for codeword in self.codewords])

    def add_parity_bit(self):
        """ Return a new Binary_code object generated by adding a parity bit to the current codewords.
        The new Binary_code will have the same number of codewords, a length higher by 1, and a minimum Hamming distance
        equal to the current one if the current one is even, or higher by one if the current one is odd."""
        new_code = Binary_code(self.length+1)
        for codeword in self.codewords:
            parity_bit = codeword.weight() % 2
            new_code.add(codeword.string()+str(parity_bit))
        assert new_code.size()==self.size()
        return new_code

    def invert(self):
        """ Return a new Binary_code object containing the bitwise inverses of all the codewords in this code."""
        new_code = Binary_code(self.length)
        for codeword in self.codewords:
            new_code.add((~codeword).string())
        return new_code

    def choose_by_bit_sum(self,low,high):
        """ Remove all codewords with bit sums outside the low-high range."""
        pass
        # TODO implement

    def reduce_by_Hamming_distance(self,low,high,min_count):
        """ Find a subset of at least min_count codewords with the Hamming distance for each pair in the low-hig range. """
        pass
        # TODO implement using clique_find.py

    # (TODO will also need to implement reduce_to_number on Binary_code, for when the code size is 1024 but we only want 700 or something - how should the codewords to use be picked?  By lower or higher weight, or whichever weight is more extreme?  Also remember to remove the all-zero codeword, or invert the code to get rid of it if we need exactly 2**k samples!)


if __name__=='__main__':
    """ If module is run directly, run tests. """
    import sys

    ### Binary_codeword tests
    if True:
        print "Testing Binary_codeword functionality..."

        # creation with different value types
        assert Binary_codeword('111').string() == '111'
        assert Binary_codeword(7).string() == '111'
        assert Binary_codeword('0b 111\n').string() == '111'
        assert Binary_codeword([True,True,True]).string() == '111'
        assert Binary_codeword(Binary_codeword('111')).string() == '111'

        # initiation length-checks/padding: 
        #   '111' has length 3, should work
        Binary_codeword('111',3,check_length=True)
        #   '111' doesn't have length 5, should raise an error
        try:                    Binary_codeword('111',5,check_length=True)
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_codeword initiation length-check didn't work!")
        #   padding a length-3 string to the same length shouldn't change anything
        assert Binary_codeword('111',3) == Binary_codeword('111')
        #   padding a length-3 string to a higher length should work
        assert Binary_codeword(7,4).string() == '0111'
        assert Binary_codeword('111',5).string() == '00111'
        #   padding a length-3 string to length 2 should raise an error
        try:                    Binary_codeword('111',2)
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_codeword initiation length-check didn't work!")

        # equality and comparison
        assert Binary_codeword('111') == Binary_codeword('111')
        assert Binary_codeword('111') != Binary_codeword('101')
        assert Binary_codeword('111') != Binary_codeword('1110')
        assert Binary_codeword('111') != Binary_codeword('0111')
        assert Binary_codeword('101') < Binary_codeword('111')
        # bitwise operations
        assert ~Binary_codeword('000') == Binary_codeword('111')
        assert Binary_codeword('110') | Binary_codeword('011') == Binary_codeword('111')
        assert Binary_codeword('110') & Binary_codeword('011') == Binary_codeword('010')
        assert Binary_codeword('110') ^ Binary_codeword('011') == Binary_codeword('101')
        # length
        assert len(Binary_codeword('111')) == 3
        assert len(Binary_codeword('00000')) == 5
        # weight
        assert Binary_codeword('111').weight() == 3
        assert Binary_codeword('001').weight() == 1
        assert Binary_codeword('00000').weight() == 0
        # string and list representations
        assert Binary_codeword('111').string() == '111'
        assert Binary_codeword('111').list() == [1,1,1]
        
        # Hamming distance
        assert Hamming_distance(Binary_codeword('000'),Binary_codeword('000')) == 0
        assert Hamming_distance(Binary_codeword('111'),Binary_codeword('000')) == 3
        assert Hamming_distance(Binary_codeword('101'),Binary_codeword('000')) == 2
        assert Hamming_distance(Binary_codeword('101'),Binary_codeword('010')) == 3
        # comparing bitstrings of different lengths should fail
        try:               Hamming_distance(Binary_codeword('000'),Binary_codeword('0000'))
        except ValueError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:              raise BinaryCodeError("Hamming distance between Binary_codewords of different lengths worked!")

        print "...DONE"

    ### Binary_code tests
    if True:
        print "Testing Binary_code functionality..."

        # basic creation from a list and all the properties
        B = Binary_code(3,['110','101','011','000'])
        assert B.length == 3
        assert B.size() == 4
        assert B.find_Hamming_distance_range() == (2,2)
        assert B.find_bit_sum_counts() == [(0,1), (2,3)]
        assert B.total_bit_sum() == 6

        # adding and removing a codeword
        B = Binary_code(3,['110','101','011','000'])
        C = Binary_code(3,B.codewords)
        C.add('111')
        assert C.length == B.length
        assert C.size() == B.size() + 1
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == B.find_bit_sum_counts() + [(3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3
        C.remove('110')
        assert C.length == B.length
        assert C.size() == B.size()
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == [(0,1), (2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3 - 2

        # removing the all-zero codeword should return 1 if it was there and 0 if it wasn't, and never fail
        assert C.remove_all_zero_codeword() == 1
        assert C.length == B.length
        assert C.size() == B.size() - 1
        assert C.find_Hamming_distance_range() == (1,2)
        assert C.find_bit_sum_counts() == [(2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3 - 2
        assert C.remove_all_zero_codeword() == 0
        assert C.length == B.length
        assert C.size() == B.size() - 1

        # inversion
        B = Binary_code(3,['110','101','011','000'])
        B_ = B.invert()
        assert B_.length == B.length
        assert B_.size() == B_.size()
        assert B_.find_Hamming_distance_range() == B.find_Hamming_distance_range()
        assert B_.find_bit_sum_counts() == sorted([(B.length-w,n) for (w,n) in B.find_bit_sum_counts()])
        assert B_.total_bit_sum() == B.size() * B.length - B.total_bit_sum()

        # adding parity bit
        D = Binary_code(2,['11','10','01','00'])
        assert D.length == 2
        assert D.size() == 4
        assert D.find_Hamming_distance_range() == (1,2)
        assert D.find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        assert D.total_bit_sum() == 4
        E = D.add_parity_bit()
        assert E.length == D.length + 1
        assert E.size() == D.size()
        assert E.find_Hamming_distance_range() == (2,2)
        assert E.find_bit_sum_counts() == [(0,1), (2,3)]
        assert E.total_bit_sum() == D.total_bit_sum() + 2

        # check that creation fails if the length is wrong
        try:                    B = Binary_code(4,['110','101','011','000'])
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code creation length-check didn't work!")
        # check that creation fails if the expected count is wrong
        try:                    B = Binary_code(3,['110','101','011','000'],expected_count=5)
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code creation expected count check didn't work!")
        # check that add/remove fails if the length is wrong
        B = Binary_code(3,['110','101','011','000'])
        try:                    B.add('1111')
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code.add(val) length-check didn't work!")
        try:                    B.remove('1111')
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code.remove(val) length-check didn't work!")
        try:                    B.remove('111')
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code.remove(val) membership-check didn't work!")

        # check that creation fails with inexistent method keyword
        try:                    B = Binary_code(4,['110','101','011','000'],method='random')
        except BinaryCodeError: pass    # this SHOULD raise an error; complain if it doesn't!
        else:                   raise BinaryCodeError("Binary_code creation length-check didn't work!")

        # check creation from matrix generator file
        infile1 = 'error-correcting_codes/19-10-5_generator'
        try:            B19 = Binary_code(19,val=infile1,method='matrixfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run matrix generator test."%infile1)
        assert B19.find_bit_sum_counts() == [(0,1), (5,30), (6,64), (7,90), (8,150), (9,180), (10,168), (11,156), (12,104), (13,46), (14,24), (15,10), (16,1)]
        B20 = B19.add_parity_bit()
        assert B20.find_bit_sum_counts() == [(0, 1), (6, 94), (8, 240), (10, 348), (12, 260), (14, 70), (16, 11)]
        # could also check Hamming distances, but those take a long time!
        #assert B19.find_Hamming_distance_range() == (5,16)
        #assert B20.find_Hamming_distance_range() == (6,16)

        # TODO check creation from code list file?  I don't have one right now - could get rid of that option, or else add a print_to_file method, use that on one of the other codes and then test reading that in.  Or just not test it - it's just reading a file, after all.

        print "...DONE"

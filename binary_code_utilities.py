#! /usr/bin/env python2
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

import sys
from collections import defaultdict
from itertools import combinations, permutations
from numpy import array, dot
import bitstring
import unittest

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
        # MAYBE-TODO in order to make this absolutely right, I should make sure the bitstrings are immutable... How?
        # MAYBE-TODO is string comparison really what I want here?  How about when the lengths are different? Should 
        #   bitstrings with different lengths even be comparable?  I suppose they should just so I can sort stuff and 
        #   get a consistent result. Possibly just sorting by length first would be better, but it doesn't matter much.
        #   As long as identity works correctly and there's SOME reproducible sorting, we're fine.
        # MAYBE-TODO should I also implement comparison to strings and such, by doing binary_codeword(other).string()?
        #   That conversion should probably be taken care of in the function performing the comparison.
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

def bit_change_count(val1,val2):
    """ Given two binary strings, return A) the count of bits that were 0 in val1 and 1 in val2, and B) the reverse. """
    count_0_to_1 = (~val1)&val2
    count_1_to_0 = val1&(~val2)
    return count_0_to_1.weight(), count_1_to_0.weight()


### Binary code (set of binary strings) representation

# MAYBE-TODO figure out a naming that won't confuse people!  (Or me!)  Mathematically a "code" is a set of codewords (or a method of encoding things), but IRL a "code" can be either a method of encoding things or just a string, so code/codeword is confusing and even I'm using them wrong!

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
        OUTFILE = open(outfile,'w')
        for codeword in self.codewords:
            OUTFILE.write(codeword.string()+'\n')
        OUTFILE.close()

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

    # MAYBE-TODO I could make one or both of these be the initialization signature instead, but who cares
    # MAYBE-TODO could add the minimum Hamming distance to this?
    def __str__(self):      return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())
    def __repr__(self):     return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())

    # Implementing eq, ne and hashing based on codeword sets; no ge/le comparison (you can't sort sets)
    # NOTE: comparison and hashing are related and need to match!  See notes in Binary_codeword.
    def __eq__(self,other):     return self.codewords == other.codewords
    def __ne__(self,other):     return self.codewords != other.codewords
    def __hash__(self):         return hash(self.codewords)

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
        """ Remove all codewords with bit sums outside the low-high range. If high is -1, don't apply an upper bound. """
        if high==-1:    high = self.find_bit_sum_counts()[-1][0]
        new_codewords = set()
        for codeword in self.codewords:
            if low <= codeword.weight() <= high:  new_codewords.add(codeword)
        self.codewords = new_codewords

    def give_N_codeword_list_by_bit_sum(self,N,remove_low=False):
        """ Sort all the codewords by bit sum (depends on remove_low), return the first N as a list.  
        If remove_low is False (default), we want to remove high-sum codewords, so sort normally and take N first ones.
        If remove_low is True, we want to remove low-sum codewords, so sort reverse and take N first (highest) ones.
        Sort by string value first to make sure the result is reproducible (so codewords with the same bit-sum are sorted
        lexicographically).  Raise an error if N is higher than current code size. """
        if N>self.size():
            raise BinaryCodeError("Cannot reduce the code to %s elements, it's already only %s!"%(N,self.size()))
        # Grab the codewords as a list; Sort normally (i.e. by string) just to make sure the result is reproducible.
        codewords = sorted(list(self.codewords))
        # Sort the codewords by weight and then take the first N (see docstring)
        codewords.sort(key = lambda codeword: codeword.weight(), reverse=remove_low)
        return codewords[:N]

    def clonality_count_conflicts(self, permitted_1_to_0_changes=0, permitted_0_to_1_changes=0, count_self_conflict=False):
        """ Simple clonality conflict count.  Return a (conflict_count: codeword_set) dictionary.
        Go over all combinations of codewords A,B,C in the code, and whenever binary_or(A,B)=C (or is close enough to C  
         according to the permitted_0_to_1_changes and permitted_1_to_0_changes arguments, where 0_to_1 means a 0 in C), 
         add one conflict count to all three codewords. 
        The count_self_conflict argument specifies whether A+B=A is considered as a problem the same way A+B=C is.
        Codewords with 0 conflicts are guaranteed not to generate problems, but nothing very useful can be said about
         how many and which of the ones with low counts could be added to that set without generating clonality issues.
        """
        # TODO this is really insanely slow!!  There's probably a way to speed it up...
        #   * can definitely speed up the version with default arguments by just going over all A,B pairs and checking for set membership of the result in the set of codewords.  DONE
        #   * when the arguments AREN'T 0, could speed it up (at the cost of memory) by generating a set of all disallowed clonality-products based on the set of all codewords, and then going over all A,B combinations and checking for membership in that set...  Yeah, seems like a better idea, doesn't it. 
        #       so TODO need to implement a way of making a set of all possible changed versions of a codeword, given the number of allowed 0-to-1 and 1-to-0 changes. All right.
        codeword_to_conflict_count = dict([(codeword,0) for codeword in self.codewords])
        if permitted_1_to_0_changes==0 and permitted_0_to_1_changes==0:
            for A,B in combinations(self.codewords,2):
                clonality_signature = A|B
                # first check the A+B=A special case:  if count_self_conflict is true, give A and B one conflict count 
                #   (not 2 and 1 as the stadard case would), otherwise do go on to the A+B=C check.
                if clonality_signature in [A,B]:
                    if count_self_conflict:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
                # now check the other codewords (note that this won't trigger for A+B=A because it's an elif)
                elif clonality_signature in self.codewords:
                    for codeword in A,B,clonality_signature:
                        codeword_to_conflict_count[codeword] += 1
        else:
            for A,B,C in permutations(self.codewords,3):
                # the order of A,B doesn't matter, so skip half; the order of C does.
                if A<B:     continue
                clonality_signature = A|B
                changes_0_to_1, changes_1_to_0 = bit_change_count(C,A|B)
                if changes_0_to_1 <= permitted_0_to_1_changes and changes_1_to_0 <= permitted_1_to_0_changes:
                    for codeword in A,B,C:
                        codeword_to_conflict_count[codeword] += 1
            # Do we want to count A+B=A as a problem, like A+B=C?  Optional separate routine.
            if count_self_conflict:
                for A,B in permutations(self.codewords,2):
                    clonality_signature = A|B
                    changes_0_to_1, changes_1_to_0 = bit_change_count(A,A|B)
                    if changes_0_to_1 <= permitted_0_to_1_changes and changes_1_to_0 <= permitted_1_to_0_changes:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
        # now generate the final conflict_count:codeword_set dictionary from the codeword:conflict_count one
        conflict_count_to_codeword_set = defaultdict(lambda: set())
        for (codeword,conflict_count) in codeword_to_conflict_count.iteritems():
            conflict_count_to_codeword_set[conflict_count].add(codeword)
        # it's not good to return defaultdicts - just transform to a dict
        return dict(conflict_count_to_codeword_set)
        # TODO I suppose I should test it on the extended Wagner code to see if my results are the same as Goodman2009...
        #   remember to remove the all-zero codeword first!

    def clonality_naive_reduction(self, permitted_1_to_0_changes=0, permitted_0_to_1_changes=0, count_self_conflict=False):
        """ Naive solution of the clonality problem, as in Goodman2009: return a set of codewords with 0 conflicts, 
         as given by clonality_count_conflicts with the same arguments (see that function's docstring for more)."""
        conflict_count_to_codeword_set = self.clonality_count_conflicts(permitted_1_to_0_changes, 
                                                                        permitted_0_to_1_changes, count_self_conflict)
        try:                return conflict_count_to_codeword_set[0]
        except KeyError:    return set()

    def reduce_by_Hamming_distance(self,low,high,min_count):
        """ Find a subset of at least min_count codewords with the Hamming distance for each pair in the low-hig range. """
        pass
        # MAYBE-TODO implement using clique_find.py

    # MAYBE-TODO implement some other options for reducing the set to one without clonality issues?
    #  * Any other sensible algorithms for doing this?  See Clonality section of notes.txt. - I had some new ideas...


class Testing_Binary_codeword(unittest.TestCase):
    """ Testing Binary_codeword functionality. """

    def test_creation_value_types(self):
        """ Test Binary_codeword instance initiation with different value types (string, int, list, Binary_codeword) """
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword(7).string(), '111')
        self.assertEqual(Binary_codeword('0b 111\n').string(), '111')
        self.assertEqual(Binary_codeword([True,True,True]).string(), '111')
        self.assertEqual(Binary_codeword(Binary_codeword('111')).string(), '111')

    def test_initiation_length_check_and_padding(self):
        """ Test if Binary_codeword initiation deals with codeword length correctly (length check and padding). """
        #   '111' has length 3, should work (no assert, just making sure it doesn't raise an exception)
        Binary_codeword('111',3,check_length=True)
        #   '111' doesn't have length 5, should raise an error
        self.assertRaises(BinaryCodeError, Binary_codeword, '111', 5, check_length=True)
        #   padding a length-3 string to the same length shouldn't change anything
        self.assertEqual(Binary_codeword('111',3), Binary_codeword('111'))
        #   padding a length-3 string to a higher length should work
        self.assertEqual(Binary_codeword(7,4).string(), '0111')
        self.assertEqual(Binary_codeword('111',5).string(), '00111')
        #   padding a length-3 string to length 2 should raise an error
        self.assertRaises(BinaryCodeError, Binary_codeword, '111', 2)

    def test_equality_and_comparison_operators(self):
        self.assertTrue(Binary_codeword('111') == Binary_codeword('111'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('101'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('1110'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('0111'))
        self.assertTrue(Binary_codeword('101') < Binary_codeword('111'))

    def test_bitwise_operations(self):
        self.assertEqual(~Binary_codeword('000'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') | Binary_codeword('011'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') & Binary_codeword('011'), Binary_codeword('010'))
        self.assertEqual(Binary_codeword('110') ^ Binary_codeword('011'), Binary_codeword('101'))
        # comparing bitstrings of different lengths should fail
        self.assertRaises(ValueError, Hamming_distance, Binary_codeword('000'), Binary_codeword('0000'))

    def test_length_calculation(self):
        self.assertEqual(len(Binary_codeword('111')), 3)
        self.assertEqual(len(Binary_codeword('00000')), 5)

    def test_weight_calculation(self):
        self.assertEqual(Binary_codeword('111').weight(), 3)
        self.assertEqual(Binary_codeword('001').weight(), 1)
        self.assertEqual(Binary_codeword('00000').weight(), 0)

    def test_list_string_representations(self):
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword('111').list(), [1,1,1])


class Testing_other_functions(unittest.TestCase):
    """ Testing functions that aren't part of either of the main classes."""
    
    def test_Hamming_distance_calculation(self):
        self.assertEqual(Hamming_distance(Binary_codeword('000'),Binary_codeword('000')), 0)
        self.assertEqual(Hamming_distance(Binary_codeword('111'),Binary_codeword('000')), 3)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('000')), 2)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('010')), 3)

    def test_bit_change_count(self):
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('000')) == (0,0)
        assert bit_change_count(Binary_codeword('111'),Binary_codeword('000')) == (0,3)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('111')) == (3,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('000')) == (0,2)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('101')) == (2,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('010')) == (1,2)
        assert bit_change_count(Binary_codeword('010'),Binary_codeword('101')) == (2,1)


class Testing_Binary_code(unittest.TestCase):
    """ Testing Binary_code functionality. """

    # MAYBE-TODO convert all the assert statements to self.assertEqual or self.assertTrue or such? 
    #   That's the way unittest functions should be written, but the current version works too...
    #   I could probably just use nosetest if I wanted - that catches normal assert statements too. 

    def test_creation_from_list_and_properties(self):
        B = Binary_code(3,['110','101','011','000'])
        assert B.length == 3
        assert B.size() == 4
        assert B.find_Hamming_distance_range() == (2,2)
        assert B.find_bit_sum_counts() == [(0,1), (2,3)]
        assert B.total_bit_sum() == 6
        # check that creation fails if the length is wrong
        self.assertRaises(BinaryCodeError, Binary_code, 4, ['110','101','011','000'])
        # check that creation fails if the expected count is wrong
        self.assertRaises(BinaryCodeError, Binary_code, 3, ['110','101','011','000'], expected_count=5)
        # check that creation fails with inexistent method keyword
        self.assertRaises(BinaryCodeError, Binary_code, 4, ['110','101','011','000'], method='random')

    def test_codeword_add_remove(self):
        B = Binary_code(3,['110','101','011','000'])
        C = Binary_code(3,B.codewords)
        # add an element, verify that the Binary_code properties are correct
        C.add('111')
        assert C.length == B.length
        assert C.size() == B.size() + 1
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == B.find_bit_sum_counts() + [(3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3
        # remove an element, verify that the Binary_code properties are correct
        C.remove('110')
        assert C.length == B.length
        assert C.size() == B.size()
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == [(0,1), (2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3 - 2
        # remove should fail if the element wasn't there
        self.assertRaises(BinaryCodeError, B.remove, '111')
        # add/remove should fail if the length is wrong
        self.assertRaises(BinaryCodeError, B.add, '1111')
        self.assertRaises(BinaryCodeError, B.remove, '1111')

    def test_remove_all_zero_codeword(self):
        B = Binary_code(3,['111','101','011','000'])
        C = Binary_code(3,B.codewords)
        assert C.remove_all_zero_codeword() == 1
        assert C.length == B.length
        assert C.size() == B.size() - 1
        assert C.find_Hamming_distance_range() == (1,2)
        assert C.find_bit_sum_counts() == [(2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum()
        # try removing the all-zero codeword again - should have no effect
        assert C.remove_all_zero_codeword() == 0
        assert C.length == B.length
        assert C.size() == B.size() - 1

    def test_codeword_inversion(self):
        B = Binary_code(3,['110','101','011','000'])
        B_ = B.invert()
        assert B_.length == B.length
        assert B_.size() == B_.size()
        assert B_.find_Hamming_distance_range() == B.find_Hamming_distance_range()
        assert B_.find_bit_sum_counts() == sorted([(B.length-w,n) for (w,n) in B.find_bit_sum_counts()])
        assert B_.total_bit_sum() == B.size() * B.length - B.total_bit_sum()

    def test_adding_parity_bit(self):
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

    def test_choose_by_bit_sum(self):
        D = Binary_code(2,['11','10','01','00'])
        assert D.length == 2
        assert D.size() == 4
        assert D.find_Hamming_distance_range() == (1,2)
        assert D.find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        assert D.total_bit_sum() == 4
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(0,2)
        assert D_ == D
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(0,0)
        assert D_.size() == 1
        assert D_.find_bit_sum_counts() == [(0,1)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(0,1)
        assert D_.size() == 3
        assert D_.find_bit_sum_counts() == [(0,1),(1,2)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(1,2)
        assert D_.size() == 3
        assert D_.find_bit_sum_counts() == [(1,2),(2,1)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(1,1)
        assert D_.size() == 2
        assert D_.find_bit_sum_counts() == [(1,2)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_by_bit_sum(0,-1)
        assert D_ == D
        E = Binary_code(2,['11','00'])
        E.choose_by_bit_sum(1,1)
        assert E.size() == 0

    def test_give_N_codeword_list_by_bit_sum(self):
        D = Binary_code(2,['11','10','01','00'])
        assert D.find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        # check that the highest (or lowest) bit-sum elements are removed first (the True/False argument is remove_low)
        #   (returns a Binary_codeword list, so make a Binary_code from it again to check bit sum counts)
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(4,False)).find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(4,True)).find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(3,False)).find_bit_sum_counts() == [(0,1), (1,2)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(3,True)).find_bit_sum_counts() == [(1,2), (2,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(2,False)).find_bit_sum_counts() == [(0,1), (1,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(2,True)).find_bit_sum_counts() == [(1,1), (2,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(1,False)).find_bit_sum_counts() == [(0,1)]
        assert Binary_code(2,D.give_N_codeword_list_by_bit_sum(1,True)).find_bit_sum_counts() == [(2,1)]
        # check the exact ordering of returned elements (sorted by bit sum, equal bit sum sorted lexicographically)
        [b01,b10,b11,b00] = [Binary_codeword(x) for x in ['01','10','11','00']]
        assert D.give_N_codeword_list_by_bit_sum(4,False) == [b00,b01,b10,b11]
        assert D.give_N_codeword_list_by_bit_sum(4,True) == [b11,b01,b10,b00]
        assert D.give_N_codeword_list_by_bit_sum(3,False) == [b00,b01,b10]
        assert D.give_N_codeword_list_by_bit_sum(3,True) == [b11,b01,b10]
        assert D.give_N_codeword_list_by_bit_sum(2,False) == [b00,b01]
        assert D.give_N_codeword_list_by_bit_sum(2,True) == [b11,b01]
        assert D.give_N_codeword_list_by_bit_sum(1,False) == [b00]
        assert D.give_N_codeword_list_by_bit_sum(1,True) == [b11]
        # check that it's impossible to get 5 elements from a 4-element binary code
        self.assertRaises(BinaryCodeError, D.give_N_codeword_list_by_bit_sum, 5, False)

    def test_creation_from_matrix_generator_file(self):
        # (also implicitly checks generation from a matrix object)
        infile1 = 'error-correcting_codes/19-10-5_generator'
        try:            B19 = Binary_code(19,val=infile1,method='matrixfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run matrix file test."%infile1)
        assert B19.find_bit_sum_counts() == [(0,1), (5,30), (6,64), (7,90), (8,150), (9,180), (10,168), (11,156), (12,104), (13,46), (14,24), (15,10), (16,1)]
        B20 = B19.add_parity_bit()
        assert B20.size() == 2**10
        assert B20.length == 20
        assert B20.find_bit_sum_counts() == [(0, 1), (6, 94), (8, 240), (10, 348), (12, 260), (14, 70), (16, 11)]

    def test_creation_from_code_list_file(self):
        B19 = Binary_code(19,val='error-correcting_codes/19-10-5_generator',method='matrixfile',expected_count=2**10)
        B20 = B19.add_parity_bit()
        infile2 = 'error-correcting_codes/20-10-6_list'
        try:            B20_new = Binary_code(20,val=infile2,method='listfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run list file test."%infile2)
        assert B20_new == B20

    def test_clonality_count_conflicts(self):
        # defining the binary codewords so I can use them to check that the results are right
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        [b001,b010,b100,b111] = [Binary_codeword(x) for x in ['001','010','100','111']]
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        B = Binary_code(3,[b110,b101,b011,b000])
        assert B.clonality_count_conflicts(0,0,False) == {0:set([b110,b101,b011,b000])}
        print B.clonality_count_conflicts(0,0,True)
        assert B.clonality_count_conflicts(0,0,True) == {1:set([b110,b101,b011]),3:set([b000])}
        C = Binary_code(3,[b110,b101,b011])
        assert C.clonality_count_conflicts(0,0,False) == {0:set([b110,b101,b011])}
        assert C.clonality_count_conflicts(0,0,True) == {0:set([b110,b101,b011])}
        assert C.clonality_count_conflicts(1,0,False) == {0:set([b110,b101,b011])}
        assert C.clonality_count_conflicts(0,1,False) == {3:set([b110,b101,b011])}
        D = Binary_code(2,[b11,b10,b01,b00])
        assert D.clonality_count_conflicts(0,0,False) == {0:set([b00]),1:set([b11,b10,b01])}
        assert D.clonality_count_conflicts(0,0,True) == {3:set([b00,b10,b01]),4:set([b11])}
        E = Binary_code(2,[b11,b00])
        assert E.clonality_count_conflicts(0,0,False) == {0:set([b11,b00])}
        assert E.clonality_count_conflicts(0,0,True) == {1:set([b11,b00])}
        F = Binary_code(3,[b001,b010,b100,b111])
        assert F.clonality_count_conflicts(0,0,False) == {0:set([b001,b010,b100,b111])}
        assert F.clonality_count_conflicts(0,0,True) == {1:set([b001,b010,b100]),3:set([b111])}
        assert F.clonality_count_conflicts(1,0,False) == {2:set([b001,b010,b100]),3:set([b111])}
        assert F.clonality_count_conflicts(0,1,False) == {0:set([b001,b010,b100,b111])}

    def test_clonality_naive_reduction(self):
        # this is pretty much strictly based on the test_clonality_count_conflicts function...
        # MAYBE-TODO there may be a way of avoiding code duplication here!
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        [b001,b010,b100,b111] = [Binary_codeword(x) for x in ['001','010','100','111']]
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        B = Binary_code(3,[b110,b101,b011,b000])
        assert B.clonality_naive_reduction(0,0,False) == set([b110,b101,b011,b000])
        assert B.clonality_naive_reduction(0,0,True) == set()
        C = Binary_code(3,[b110,b101,b011])
        assert C.clonality_naive_reduction(0,0,False) == set([b110,b101,b011])
        assert C.clonality_naive_reduction(0,0,True) == set([b110,b101,b011])
        assert C.clonality_naive_reduction(1,0,False) == set([b110,b101,b011])
        assert C.clonality_naive_reduction(0,1,False) == set()
        D = Binary_code(2,[b11,b10,b01,b00])
        assert D.clonality_naive_reduction(0,0,False) == set([b00])
        assert D.clonality_naive_reduction(0,0,True) == set()
        E = Binary_code(2,[b11,b00])
        assert E.clonality_naive_reduction(0,0,False) == set([b11,b00])
        assert E.clonality_naive_reduction(0,0,True) == set()
        F = Binary_code(3,[b001,b010,b100,b111])
        assert F.clonality_naive_reduction(0,0,False) == set([b001,b010,b100,b111])
        assert F.clonality_naive_reduction(0,0,True) == set()
        assert F.clonality_naive_reduction(1,0,False) == set()
        assert F.clonality_naive_reduction(0,1,False) == set([b001,b010,b100,b111])


if __name__=='__main__':
    """ If module is ran directly, run tests. """

    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()

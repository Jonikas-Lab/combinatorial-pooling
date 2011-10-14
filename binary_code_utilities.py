#! /usr/bin/env python2
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

import sys
from collections import defaultdict
import itertools
from numpy import array, dot
import bitstring
from general_utilities import invert_dict_tolists, invert_listdict_tolists
import unittest

class BinaryCodeError(Exception):
    """ Exceptions in the binary_code_utilities module."""
    pass

######### Binary string representations
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


######### basic functions on multiple Binary_codeword objects: distance, mutations, etc

def Hamming_distance(val1,val2):
    """ Given two binary strings, return the number of bits by which their binary representations differ. """
    bitwise_xor = val1^val2
    return bitwise_xor.weight()


def bit_change_count(val1,val2):
    """ Given two binary strings, return A) the count of bits that were 0 in val1 and 1 in val2, and B) the reverse. """
    count_0_to_1 = (~val1)&val2
    count_1_to_0 = val1&(~val2)
    return count_0_to_1.weight(), count_1_to_0.weight()


def _change_all_position_combinations(val, max_N_changes, allowed_change_positions=None):
    """ Just a help function for expand_by_all_mutations.
    Return a set of all new values generated by inverting (0/1) at most max_N_changes of the bits 
     of the original val argument, with changes only allowed in allowed_change_positions (all by default).
    val must be a list of 0/1 ints. """
    try:                    val[0]
    except TypeError:       raise ValueError("val has to be a list of 0/1 ints!")
    if not val[0] in [0,1]: raise ValueError("val has to be a list of 0/1 ints!")
    # if no allowed_change_positions was given, assume all positions are allowed
    if allowed_change_positions is None:
        allowed_change_positions = range(len(val))
    new_value_set = set()
    for N_changes in range(max_N_changes+1):
        if N_changes>len(allowed_change_positions): break   # this is necessary in 2.6.1 but not 2.6.5...
        for change_positions in itertools.combinations(allowed_change_positions,N_changes):
            new_val = list(val)
            for pos in change_positions:
                new_val[pos] = int(not val[pos])
            new_value_set.add(tuple(new_val))   # need to convert to tuple because lists aren't hashable
    return new_value_set


def expand_by_all_mutations(codeword_set, N_permitted_changes):
    """ Return set of all codewords that are too close to codeword_set (distant by at most the given number of changes).
    N_permitted_changes can be given either as a single value, or as a (N_1_to_0_changes, N_0_to_1_changes) tuple."""
    codewords_as_lists = [codeword.list() for codeword in codeword_set]
    expanded_codeword_set = set()

    if isinstance(N_permitted_changes,int):
        for codeword in codewords_as_lists:
            # change all possible combinations of positions (0s or 1s, all positions are allowed)
            new_val_full_set = _change_all_position_combinations(codeword,N_permitted_changes)
            expanded_codeword_set.update([Binary_codeword(new_val) for new_val in new_val_full_set])

    elif isinstance(N_permitted_changes,tuple) and len(N_permitted_changes)==2:
        (permitted_1_to_0_changes,permitted_0_to_1_changes) = N_permitted_changes
        for codeword in codewords_as_lists:
            # change all possible combinations of 1s to 0s up to the allowed limit
            list_of_1_positions = [x for x in range(len(codeword)) if codeword[x]==1]
            new_val_halfway_set = _change_all_position_combinations(codeword,permitted_1_to_0_changes,list_of_1_positions)
            # now take each of the resulting values and change all possible combinations of 0s up to the allowed limit
            for new_val in new_val_halfway_set:
                list_of_0_positions = [x for x in range(len(new_val)) if new_val[x]==0]
                new_val_full_set = _change_all_position_combinations(new_val, permitted_0_to_1_changes,list_of_0_positions)
                expanded_codeword_set.update([Binary_codeword(new_val) for new_val in new_val_full_set])
    else:
        raise ValueError("N_permitted_changes must be an int or a tuple of two ints!")

    return expanded_codeword_set


def expand_by_all_mutations_dict(codeword_set, N_permitted_changes):
    """ Return a dictionary with the keys being all the codeword too close to codeword_set, 
     and the values being the list of base codewords that that particular result codeword was too close to.
    (See expand_by_all_mutations docstring for details on what "too close" means.) """
    codeword_to_expanded_set = dict()
    for C in codeword_set:
        codeword_to_expanded_set[C] = expand_by_all_mutations([C], N_permitted_changes)
    return invert_listdict_tolists(codeword_to_expanded_set)


######### Binary code (set of binary strings) representation

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
        # MAYBE-TODO or should I make this return a new code instead? Do I want the codes to be immutable or not?
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
    # MAYBE-TODO this hashing solution is dangerous, since technically codeword sets ARE mutable... Hmmmm...
    #def __hash__(self):         return hash(frozenset(self.codewords))      

    def find_Hamming_distance_range(self):
        """ Return a tuple containing the lowest and highest Hamming distance between all codeword pairs."""
        # set start lov/high values to impossible extremes to make sure they get overwritten
        if len(self.codewords)==0:
            return None,None
        low = self.length+1
        high = 0
        for x,y in itertools.combinations(self.codewords,2):
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

    def clonality_count_conflicts(self, N_allowed_changes=(0,0), count_self_conflict=False, remove_all_zero_codeword=False,
                                  print_conflict_details=False, return_conflict_details=False, quiet=False):
        """ Simple clonality conflict count.  Return a (conflict_count: codeword_set) dictionary.
        Go over all combinations of codewords A,B,C in the code, and whenever the clonality sum A+B is close enough to C  
         according to N_allowed changes (which can be either a single number or a (1_to_0_changes, 0_to_1_changes) tuple)
         add one conflict count to all three codewords. 
        The count_self_conflict argument specifies whether A+B=A is considered as a problem the same way A+B=C is; 
         if count_self_conflict is True, you may want to also set remove_all_zero_codeword to True, 
         otherwise it'll conflict with everything (a warning will be printed).
        Codewords with 0 conflicts are guaranteed not to generate problems, but nothing very useful can be said about
         how many and which of the ones with low counts could be added to that set without generating clonality issues.
        The last two arguments allow printing and/or returning of detailed conflict data - a set of (A,B,A|B,C,s,n) 
         tuples, where A and B are two codewords in the code, A|B is their clonality result, C is the set of codewords
          that A|B conflicts with (i.e. is too close to, based on N_allowed_changes), s is 'self-conflict' if it was 
          a self-conflict and '' otherwise, and n is N_allowed_changes (since I don't have the real N_changes available).
        """

        # deal with all-zero codeword: remove if requested, print warning if it's False and count_self_conflict is True
        if remove_all_zero_codeword:    self.remove_all_zero_codeword()
        elif count_self_conflict and (Binary_codeword('0'*self.length,length=self.length) in self.codewords) and not quiet:
            print("Warning: you're running a clonality conflict check with count_self_conflict turned on, and your code "
                  +"contains the all-zero codeword - be aware that it's going to generate a clonality conflict with "
                  +"EVERYTHING. Set the remove_all_zero_codeword argument to True if you'd like to prevent that; "
                  +"set the quiet argument to True to silence this message.")

        # set up conflict-count dictionary, with a 0 for each codeword
        codeword_to_conflict_count = dict([(codeword,0) for codeword in self.codewords])
        if return_conflict_details:     all_conflict_details = set()

        # TODO implement print_conflict_details and return_conflict_details!

        ### Special case just for 0 allowed_changes, because the normal way is SLOW
        # MAYBE-TODO it would be better code if the special case wasn't here... But it is faster than the general case.
        # MAYBE-TODO add option to force using the general case even for 0 changes, to make totally sure it's the same?
        if N_allowed_changes in [0, (0,0)]:
            for A,B in itertools.combinations(self.codewords,2):
                clonality_result = A|B
                conflict_details = None
                # if clonality_result is A or B, what to do depends on if count_self_conflict is True
                if clonality_result in [A,B]:
                    if count_self_conflict:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
                        conflict_details = (A,B,clonality_result,'self-conflict',N_allowed_changes)
                # if clonality_result is an existing codeword (that's not A or B), add a conflict count
                elif clonality_result in self.codewords: 
                    for codeword in [A,B,clonality_result]:
                        codeword_to_conflict_count[codeword] += 1
                    conflict_details = (A,B,clonality_result,clonality_result,'',N_allowed_changes)
                if conflict_details and print_conflict_details:     print conflict_details
                if conflict_details and return_conflict_details:    all_conflict_details.add(conflict_details)

        ### Standard case, for when the allowed change count arguments aren't 0: 
        #   pre-calculate a pool of all illegal clonality results based on all the codewords and check against that.
        # MAYBE-TODO why is this still so much slower than the special case version, even with 0,0 arguments?  
        #   Is it the set operations? Ir is it that much slower at all, really?... Do I care enough to fix it?
        else:
            expanded_conflict_values = expand_by_all_mutations_dict(self.codewords, N_allowed_changes)
            for A,B in itertools.combinations(self.codewords,2):
                clonality_result = A|B
                conflict_details = None
                if clonality_result in expanded_conflict_values: 
                    # if the list of base codewords for the conflict contains codewords other than A and B, 
                    #   register a conflict for A, B and all the base codewords (doing a set union 
                    #    ensures that even if A/B were among the base codewords, they only get one conflict)
                    if len(expanded_conflict_values[clonality_result]-set([A,B]))>0:
                        for codeword in expanded_conflict_values[clonality_result].union(set([A,B])):
                            codeword_to_conflict_count[codeword] += 1
                        conflict_set = expanded_conflict_values[clonality_result] - set([A,B])
                        conflict_details = (A,B,clonality_result,conflict_set,'',N_allowed_changes)
                    # if the list of base codewords for the conflict was only A/B, 
                    #   only register a conflict is count_self_conflict is True
                    elif count_self_conflict:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
                        conflict_set = expanded_conflict_values[clonality_result] & set([A,B])
                        conflict_details = (A,B,clonality_result,conflict_set,'self-conflict',N_allowed_changes)
                if conflict_details and print_conflict_details:     print conflict_details
                if conflict_details and return_conflict_details:    all_conflict_details.add(conflict_details)

        ### generate the final conflict_count:codeword_set dictionary from the codeword:conflict_count one
        conflict_count_to_codeword_set = invert_dict_tolists(codeword_to_conflict_count)
        if return_conflict_details:     return conflict_count_to_codeword_set, all_conflict_details
        else:                           return conflict_count_to_codeword_set
        # TODO this needs more testing - some results I get are odd!  
        #   (see experiments/generating_library/1110_clonality_check_troubleshooting/notes.txt for details/TODOs/etc)


    def clonality_conflict_check(self, N_allowed_changes=(0,0), count_self_conflict=False, remove_all_zero_codeword=False,
                                  print_conflict_details=False, quiet=False):
        """ Return False if the code contains no clonality conflicts based on arguments, True otherwise.
        Passes all its arguments to clonality_count_conflicts - see that function's docstring for details."""
        # MAYBE-TODO could rewrite this to be much faster by writing it separately and having it just go on until
        #   the first conflict and then return True, instead of using clonality_count_conflicts and thus needing to
        #   go through all the combinations, but I'm not sure it's worth it.
        conflict_count_to_codeword_set = self.clonality_count_conflicts(N_allowed_changes, count_self_conflict, remove_all_zero_codeword, print_conflict_details, return_conflict_details=False, quiet=quiet)
        if conflict_count_to_codeword_set.keys() == [0]:    return False
        else:                                               return True

    def clonality_naive_reduction(self, N_allowed_changes=(0,0), count_self_conflict=False, remove_all_zero_codeword=False,
                                  print_conflict_details=False, quiet=False):
        """ Naive solution of the clonality problem, as in Goodman2009: return a set of codewords with 0 conflicts, 
         as given by clonality_count_conflicts with the same arguments (see that function's docstring for more)."""
        conflict_count_to_codeword_set = self.clonality_count_conflicts(N_allowed_changes, count_self_conflict, remove_all_zero_codeword, print_conflict_details, return_conflict_details=False, quiet=quiet)
        try:                return conflict_count_to_codeword_set[0]
        except KeyError:    return set()

    def reduce_by_Hamming_distance(self,low,high,min_count):
        """ Find a subset of at least min_count codewords with the Hamming distance for each pair in the low-high range."""
        pass
        # MAYBE-TODO implement using clique_find.py

    # MAYBE-TODO implement some other options for reducing the set to one without clonality issues?
    #  * Any other sensible algorithms for doing this?  See Clonality section of notes.txt. - I had some new ideas...


class Testing__Binary_codeword(unittest.TestCase):
    """ Testing Binary_codeword functionality. """

    def test__creation_value_types(self):
        """ Test Binary_codeword instance initiation with different value types (string, int, list, Binary_codeword) """
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword(7).string(), '111')
        self.assertEqual(Binary_codeword('0b 111\n').string(), '111')
        self.assertEqual(Binary_codeword([True,True,True]).string(), '111')
        self.assertEqual(Binary_codeword(Binary_codeword('111')).string(), '111')

    def test__initiation_length_check_and_padding(self):
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

    def test__equality_and_comparison_operators(self):
        self.assertTrue(Binary_codeword('111') == Binary_codeword('111'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('101'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('1110'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('0111'))
        self.assertTrue(Binary_codeword('101') < Binary_codeword('111'))

    def test__bitwise_operations(self):
        self.assertEqual(~Binary_codeword('000'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') | Binary_codeword('011'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') & Binary_codeword('011'), Binary_codeword('010'))
        self.assertEqual(Binary_codeword('110') ^ Binary_codeword('011'), Binary_codeword('101'))
        # comparing bitstrings of different lengths should fail
        self.assertRaises(ValueError, Hamming_distance, Binary_codeword('000'), Binary_codeword('0000'))

    def test__length_calculation(self):
        self.assertEqual(len(Binary_codeword('111')), 3)
        self.assertEqual(len(Binary_codeword('00000')), 5)

    def test__weight_calculation(self):
        self.assertEqual(Binary_codeword('111').weight(), 3)
        self.assertEqual(Binary_codeword('001').weight(), 1)
        self.assertEqual(Binary_codeword('00000').weight(), 0)

    def test__list_string_representations(self):
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword('111').list(), [1,1,1])


class Testing__other_functions(unittest.TestCase):
    """ Testing functions that aren't part of either of the main classes."""
    
    def test__Hamming_distance_calculation(self):
        self.assertEqual(Hamming_distance(Binary_codeword('000'),Binary_codeword('000')), 0)
        self.assertEqual(Hamming_distance(Binary_codeword('111'),Binary_codeword('000')), 3)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('000')), 2)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('010')), 3)

    def test__bit_change_count(self):
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('000')) == (0,0)
        assert bit_change_count(Binary_codeword('111'),Binary_codeword('000')) == (0,3)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('111')) == (3,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('000')) == (0,2)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('101')) == (2,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('010')) == (1,2)
        assert bit_change_count(Binary_codeword('010'),Binary_codeword('101')) == (2,1)

    def test__change_all_position_combinations(self):
        assert _change_all_position_combinations([1,1], 0, None) == set([(1,1)])
        assert _change_all_position_combinations([1,1], 1, None) == set([(1,1),(0,1),(1,0)])
        assert _change_all_position_combinations([1,1], 2, None) == set([(1,1),(0,1),(1,0),(0,0)])
        assert _change_all_position_combinations([1,1], 10, None) == set([(1,1),(0,1),(1,0),(0,0)])
        assert _change_all_position_combinations([1,1], 1, [0]) == set([(1,1),(0,1)])
        assert _change_all_position_combinations([1,1], 1, [1]) == set([(1,1),(1,0)])
        assert _change_all_position_combinations([1,1], 1, [0,1]) == set([(1,1),(0,1),(1,0)])
        assert _change_all_position_combinations([0,0], 0, None) == set([(0,0)])
        assert _change_all_position_combinations([0,0], 1, None) == set([(0,0),(0,1),(1,0)])

    def test__expand_by_all_mutations(self):
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        all_test_codewords = set([b11,b10,b01,b00])
        single_value_results = dict()
        ### setup of values
        # N_permitted changes as a single number
        single_value_results[(b11,0)] = set([b11])
        single_value_results[(b11,1)] = set([b11,b10,b01])
        single_value_results[(b11,2)] = set([b11,b10,b01,b00])
        single_value_results[(b00,0)] = set([b00])
        single_value_results[(b00,1)] = set([b00,b10,b01])
        single_value_results[(b00,2)] = set([b00,b10,b01,b11])
        single_value_results[(b01,0)] = set([b01])
        single_value_results[(b01,1)] = set([b01,b00,b11])
        single_value_results[(b01,2)] = set([b01,b00,b11,b10])
        single_value_results[(b10,0)] = set([b10])
        single_value_results[(b10,1)] = set([b10,b00,b11])
        single_value_results[(b10,2)] = set([b10,b00,b11,b01])
        # N_permitted changes as a (N_permitted_1_to_0,N_permitted_0_to_1) tuple
        for N_changes in [(0,0),(0,1),(0,2)]:   single_value_results[(b11,N_changes)] = set([b11])
        for N_changes in [(1,0),(1,1),(1,2)]:   single_value_results[(b11,N_changes)] = set([b11,b01,b10])
        for N_changes in [(2,0),(2,1),(2,2)]:   single_value_results[(b11,N_changes)] = set([b11,b01,b10,b00])
        for N_changes in [(0,0),(1,0),(2,0)]:       single_value_results[(b00,N_changes)] = set([b00])
        for N_changes in [(0,1),(1,1),(2,1)]:       single_value_results[(b00,N_changes)] = set([b00,b01,b10])
        for N_changes in [(0,2),(1,2),(2,2)]:       single_value_results[(b00,N_changes)] = set([b00,b01,b10,b11])
        for N_changes in [(0,0)]:                       single_value_results[(b01,N_changes)] = set([b01])
        for N_changes in [(0,1),(0,2)]:                 single_value_results[(b01,N_changes)] = set([b01,b11])
        for N_changes in [(1,0),(2,0)]:                 single_value_results[(b01,N_changes)] = set([b01,b00])
        for N_changes in [(1,1),(1,2),(2,1),(2,2)]:     single_value_results[(b01,N_changes)] = set([b01,b00,b11,b10])
        for N_changes in [(0,0)]:                           single_value_results[(b10,N_changes)] = set([b10])
        for N_changes in [(0,1),(0,2)]:                     single_value_results[(b10,N_changes)] = set([b10,b11])
        for N_changes in [(1,0),(2,0)]:                     single_value_results[(b10,N_changes)] = set([b10,b00])
        for N_changes in [(1,1),(1,2),(2,1),(2,2)]:         single_value_results[(b10,N_changes)] = set([b10,b00,b11,b01])
        # now do the actual testing based on the defined values
        N_changes_values = [0,1,2]
        for N_changes in N_changes_values+list(itertools.product(N_changes_values,repeat=2)):
            for test_set_size in range(5):
                for test_codewords in itertools.combinations(all_test_codewords,test_set_size):
                    # test expand_by_all_mutations - simply a union of the result sets
                    all_result_sets = [single_value_results[(c,N_changes)] for c in test_codewords]
                    assert expand_by_all_mutations(test_codewords, N_changes) == set().union(*all_result_sets)
                    # text expand_by_all_mutations_dict - an inverted dictionary of the above, pretty much
                    all_result_dict = dict([(c,single_value_results[(c,N_changes)]) for c in test_codewords])
                    result_to_base_dict = invert_listdict_tolists(all_result_dict)
                    assert expand_by_all_mutations_dict(test_codewords, N_changes) == result_to_base_dict




class Testing__Binary_code(unittest.TestCase):
    """ Testing Binary_code functionality. """

    # MAYBE-TODO convert all the assert statements to self.assertEqual or self.assertTrue or such? 
    #   That's the way unittest functions should be written, but the current version works too...
    #   I could probably just use nosetest if I wanted - that catches normal assert statements too. 

    def test__creation_from_list_and_properties(self):
        for l in [0,1,5,100]:
            B = Binary_code(l,[])
            assert B.length == l
            assert B.size() == 0
            assert B.find_Hamming_distance_range() == (None,None)
            assert B.find_bit_sum_counts() == []
            assert B.total_bit_sum() == 0
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

    def test__codeword_add_remove(self):
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

    def test__remove_all_zero_codeword(self):
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

    def test__codeword_inversion(self):
        B = Binary_code(3,['110','101','011','000'])
        B_ = B.invert()
        assert B_.length == B.length
        assert B_.size() == B_.size()
        assert B_.find_Hamming_distance_range() == B.find_Hamming_distance_range()
        assert B_.find_bit_sum_counts() == sorted([(B.length-w,n) for (w,n) in B.find_bit_sum_counts()])
        assert B_.total_bit_sum() == B.size() * B.length - B.total_bit_sum()

    def test__adding_parity_bit(self):
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

    def test__choose_by_bit_sum(self):
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

    def test__give_N_codeword_list_by_bit_sum(self):
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

    def test__creation_from_matrix_generator_file(self):
        # (also implicitly checks generation from a matrix object)
        infile1 = 'error-correcting_codes/19-10-5_generator'
        try:            B19 = Binary_code(19,val=infile1,method='matrixfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run matrix file test."%infile1)
        assert B19.find_bit_sum_counts() == [(0,1), (5,30), (6,64), (7,90), (8,150), (9,180), (10,168), (11,156), (12,104), (13,46), (14,24), (15,10), (16,1)]
        B20 = B19.add_parity_bit()
        assert B20.size() == 2**10
        assert B20.length == 20
        assert B20.find_bit_sum_counts() == [(0, 1), (6, 94), (8, 240), (10, 348), (12, 260), (14, 70), (16, 11)]

    def test__creation_from_code_list_file(self):
        B19 = Binary_code(19,val='error-correcting_codes/19-10-5_generator',method='matrixfile',expected_count=2**10)
        B20 = B19.add_parity_bit()
        infile2 = 'error-correcting_codes/20-10-6_list'
        try:            B20_new = Binary_code(20,val=infile2,method='listfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run list file test."%infile2)
        assert B20_new == B20

    def test__clonality_count_conflicts_and_derivatives(self):
        """ Tests clonality_count_conflicts, clonality_naive_reduction, and clonality_conflict_check (since those 
        functions are trivially derived from clonality_count_conflicts and are easiest to test with the same setup)"""
        # defining the binary codewords so I can use them to check that the results are right
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        [b001,b010,b100,b111] = [Binary_codeword(x) for x in ['001','010','100','111']]
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        data_and_outputs = []
        B = Binary_code(3,[b110,b101,b011,b000])
        data_and_outputs.append((B,0,False, {0:set([b110,b101,b011,b000])} ))
        data_and_outputs.append((B,0,True,  {1:set([b110,b101,b011]),3:set([b000])} ))
        data_and_outputs.append((B,(0,0),False, {0:set([b110,b101,b011,b000])} ))
        data_and_outputs.append((B,(0,0),True,  {1:set([b110,b101,b011]),3:set([b000])} ))
        C = Binary_code(3,[b110,b101,b011])
        data_and_outputs.append((C,0,False, {0:set([b110,b101,b011])} ))
        data_and_outputs.append((C,0,True,  {0:set([b110,b101,b011])} ))
        data_and_outputs.append((C,1,False, {3:set([b110,b101,b011])} ))
        data_and_outputs.append((C,(0,0),False, {0:set([b110,b101,b011])} ))
        data_and_outputs.append((C,(0,0),True,  {0:set([b110,b101,b011])} ))
        data_and_outputs.append((C,(1,0),False, {0:set([b110,b101,b011])} ))
        data_and_outputs.append((C,(0,1),False, {3:set([b110,b101,b011])} ))
        D = Binary_code(2,[b11,b10,b01,b00])
        data_and_outputs.append((D,0,False, {0:set([b00]),1:set([b11,b10,b01])} ))
        data_and_outputs.append((D,0,True,  {3:set([b00,b10,b01]),4:set([b11])} ))
        data_and_outputs.append((D,(0,0),False, {0:set([b00]),1:set([b11,b10,b01])} ))
        data_and_outputs.append((D,(0,0),True,  {3:set([b00,b10,b01]),4:set([b11])} ))
        E = Binary_code(2,[b11,b00])
        data_and_outputs.append((E,0,False, {0:set([b11,b00])} ))
        data_and_outputs.append((E,0,True,  {1:set([b11,b00])} ))
        data_and_outputs.append((E,(0,0),False, {0:set([b11,b00])} ))
        data_and_outputs.append((E,(0,0),True,  {1:set([b11,b00])} ))
        F = Binary_code(3,[b001,b010,b100,b111])
        data_and_outputs.append((F,0,False, {0:set([b001,b010,b100,b111])} ))
        data_and_outputs.append((F,0,True,  {1:set([b001,b010,b100]),3:set([b111])} ))
        data_and_outputs.append((F,1,False, {2:set([b001,b010,b100]),3:set([b111])} ))
        data_and_outputs.append((F,(0,0),False, {0:set([b001,b010,b100,b111])} ))
        data_and_outputs.append((F,(0,0),True,  {1:set([b001,b010,b100]),3:set([b111])} ))
        data_and_outputs.append((F,(1,0),False, {2:set([b001,b010,b100]),3:set([b111])} ))
        data_and_outputs.append((F,(0,1),False, {0:set([b001,b010,b100,b111])} ))
        for (code,N_changes,self_conflict,result) in data_and_outputs:
            assert code.clonality_count_conflicts(N_changes,self_conflict,quiet=True) == result
            if 0 in result:     expected_outcome_1 = result[0]
            else:               expected_outcome_1 = set()
            assert code.clonality_naive_reduction(N_changes,self_conflict,quiet=True) == expected_outcome_1
            expected_outcome_2 = (False if result.keys()==[0] else True)
            assert code.clonality_conflict_check(N_changes,self_conflict,quiet=True) == expected_outcome_2 
        # TODO rewrite the tests using the return_conflict_details=True option, for more clarity!


if __name__=='__main__':
    """ If module is ran directly, run tests. """

    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()

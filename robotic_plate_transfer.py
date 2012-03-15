#!/usr/bin/env python2
"""
Generate Biomek transfer files to combine N samples into M pools, using a provided binary code.

Given a number of samples, a number of pools, a binary code, and information about the physical plates used for the samples and pools, assign one codeword from the code to each sample.  Write a varying number of Biomek robot transfer files (depends on -m/-o and -M options) that when given the appropriate plates will execute the combinatorial pooling.  Also write a general output file containing all the script information, the plate/well positions and codewords for each sample, and the plate/well positions for each pool. 

Assumes the samples will be in sequential order over the provided number of same-sized input plates: [plate1 well A1, plate1 well A2, ..., plate1 well B2, .., plate1 well A1, ...].  The pools will be similarly distributed over the provided number of output plates.  The number and size of output plates does not need to have any particular relationship to the number and size of input plates.

Arbitrarily assigns codewords from the provided binary code to each sample; the pools correspond to each bit of the codeword.  Whenever a bit of the codeword for the sample is 1, that sample should be added to the corresponding pool, so a line of the form "sample_plate_name,sample_plate_well,pool_plate_name,pool_plate_well,volume" is added to the Biomek command list. 

For more details on inputs and available options, run script with -h, or see help strings in define_option_parser function.

 -- Weronika Patena, Jonikas lab, Carnegie Institution, July 2011

USAGE:  robotic_plate_transfer.py [options] outfile_base_name
        robotic_plate_transfer.py [-h] [-t] [-T]
"""

# standard libraries
import sys, os, shutil
import unittest
from collections import defaultdict
from math import ceil
from string import ascii_uppercase as letters    # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# my modules
import binary_code_utilities
from general_utilities import invert_list_to_dict, save_line_list_as_file, write_header_data

class PlateTransferError(Exception):
    """ Exception for this file (does nothing interesting)."""
    pass


class Plate_type:
    """ A plate type (384-well, 96-well, 24-well or 6-well)."""

    plate_types_rows_columns = {6: (2,3), 24: (4,6), 96: (8,12), 384: (16,24)}
    plate_sizes = sorted(plate_types_rows_columns.keys())

    def __init__(self,size):
        """ Go from the size to row and column numbers, call make_well_ID_list_and_dict."""
        self.size = size
        try:                
            self.rows, self.columns = self.__class__.plate_types_rows_columns[size]
        except KeyError:    
            raise PlateTransferError("Plate size must be in integer in %s!"%self.__class__.plate_sizes)
        assert self.rows*self.columns == self.size
        self.make_well_ID_list_and_dict()
        # TODO add an option of explicily specifying the well list (or dict) instead of just size!
        # TODO use that to add a 6-well plate that's masquarading as a 96-well plate.
    
    def make_well_ID_list_and_dict(self):
        """ Generate a well ID list so list[number]==ID, and a well ID dict so dict[ID]==number.
        The list is of the form ['A1','A2', ...] and the dictionary {'A1':1, 'A2':2, ...}."""
        self.well_ID_list = []
        for row in range(self.rows):
            self.well_ID_list.extend(['%s%d'%(letters[row],col+1) for col in range(self.columns)])
        assert len(self.well_ID_list) == self.size
        self.well_ID_dict = invert_list_to_dict(self.well_ID_list)

    def get_well_ID_from_number(self,number):
        """ Given a 0-based well number (4), return the ID (B2 for 6-well plate, A5 for 96-well plate)."""
        try:                return self.well_ID_list[number]
        except IndexError:  raise PlateTransferError("Can't get well %s from a %s-well plate!"%(number,self.size))

    def get_well_number_from_ID(self,ID):
        """ Given a well ID (B1), return the sequential 0-based well number (3 for 6-well plate, 12 for 96-well plate)."""
        try:                return self.well_ID_dict[ID]
        except IndexError:  raise PlateTransferError("Can't get well %s from a %s-well plate!"%(ID,self.size))

# initialize the defined plate types
plate_types = dict([(n,Plate_type(n)) for n in Plate_type.plate_sizes])


### Unit-tests for the classes/setup above and "general functions" below
# MAYBE-TODO do I want separate test classes for each function, or should I merge them? What's the standard? Unittest output doesn't give the class name/docstring when it fails, only the function name/docstring...

class Testing__Plate_type(unittest.TestCase):
    """ Unit-tests for the Plate_type class and methods. """
    # (all the plate types have been generated already on module import, so just test them)

    def test__get_well_ID_from_number__known_wells(self):
        # try the first, second, first-in-second-row, and last well of each plate size
        assert plate_types[6].get_well_ID_from_number(0) == 'A1'
        assert plate_types[6].get_well_ID_from_number(1) == 'A2'
        assert plate_types[6].get_well_ID_from_number(3) == 'B1'
        assert plate_types[6].get_well_ID_from_number(5) == 'B3'
        assert plate_types[24].get_well_ID_from_number(0) == 'A1'
        assert plate_types[24].get_well_ID_from_number(1) == 'A2'
        assert plate_types[24].get_well_ID_from_number(6) == 'B1'
        assert plate_types[24].get_well_ID_from_number(23) == 'D6'
        assert plate_types[96].get_well_ID_from_number(0) == 'A1'
        assert plate_types[96].get_well_ID_from_number(1) == 'A2'
        assert plate_types[96].get_well_ID_from_number(12) == 'B1'
        assert plate_types[96].get_well_ID_from_number(95) == 'H12'
        assert plate_types[384].get_well_ID_from_number(0) == 'A1'
        assert plate_types[384].get_well_ID_from_number(1) == 'A2'
        assert plate_types[384].get_well_ID_from_number(24) == 'B1'
        assert plate_types[384].get_well_ID_from_number(383) == 'P24'

    def test__get_well_ID_from_number__last_wells(self):
        # you can also use a negative list index to get wells from the end, why not
        assert plate_types[6].get_well_ID_from_number(-1) == 'B3'
        assert plate_types[24].get_well_ID_from_number(-1) == 'D6'
        assert plate_types[96].get_well_ID_from_number(-1) == 'H12'
        assert plate_types[384].get_well_ID_from_number(-1) == 'P24'

    def test__get_well_ID_from_number__bad_numbers(self):
        # basic tests with obviously wrong values
        self.assertRaises(PlateTransferError, plate_types[6].get_well_ID_from_number, 100)
        self.assertRaises(TypeError, plate_types[6].get_well_ID_from_number, 0.5)
        self.assertRaises(TypeError, plate_types[6].get_well_ID_from_number, 'A')
        self.assertRaises(TypeError, plate_types[6].get_well_ID_from_number, [1])
        # well numbers are 0-based, so there should be no well N in an N-well plate (the wells are 0..N-1)
        for plate_size in Plate_type.plate_sizes:
            self.assertRaises(PlateTransferError, plate_types[plate_size].get_well_ID_from_number, plate_size)

    def test__get_well_number_from_ID__known_wells(self):
        # try the first, second, first-in-second-row, and last well of each plate size
        assert plate_types[6].get_well_number_from_ID('A1') == 0
        assert plate_types[6].get_well_number_from_ID('A2') == 1
        assert plate_types[6].get_well_number_from_ID('B1') == 3
        assert plate_types[6].get_well_number_from_ID('B3') == 5
        assert plate_types[24].get_well_number_from_ID('A1') == 0
        assert plate_types[24].get_well_number_from_ID('A2') == 1
        assert plate_types[24].get_well_number_from_ID('B1') == 6
        assert plate_types[24].get_well_number_from_ID('D6') == 23
        assert plate_types[96].get_well_number_from_ID('A1') == 0
        assert plate_types[96].get_well_number_from_ID('A2') == 1
        assert plate_types[96].get_well_number_from_ID('B1') == 12
        assert plate_types[96].get_well_number_from_ID('H12') == 95
        assert plate_types[384].get_well_number_from_ID('A1') == 0
        assert plate_types[384].get_well_number_from_ID('A2') == 1
        assert plate_types[384].get_well_number_from_ID('B1') == 24
        assert plate_types[384].get_well_number_from_ID('P24') == 383
        
    def test__get_well_number__fromID__bad_numbers(self):
        """ For each size, test a row and column that shouldn't exist). """
        # MAYBE-TODO make.get_well_number_from_ID check for this explicitly, raise PlateTransferError instead of KeyError?
        self.assertRaises(KeyError, plate_types[6].get_well_number_from_ID, 'A4')
        self.assertRaises(KeyError, plate_types[6].get_well_number_from_ID, 'C1')
        self.assertRaises(KeyError, plate_types[24].get_well_number_from_ID, 'A7')
        self.assertRaises(KeyError, plate_types[24].get_well_number_from_ID, 'E1')
        self.assertRaises(KeyError, plate_types[96].get_well_number_from_ID, 'A13')
        self.assertRaises(KeyError, plate_types[96].get_well_number_from_ID, 'I1')
        self.assertRaises(KeyError, plate_types[384].get_well_number_from_ID, 'A25')
        self.assertRaises(KeyError, plate_types[384].get_well_number_from_ID, 'Q1')

    def test__creating_bad_plate_types(self):
        # Shouldn't ever be able to create a 0-well plate
        self.assertRaises(PlateTransferError,Plate_type,0)
        # Shouldn't be able to create a 10-well plate, it wasn't defined
        self.assertRaises(PlateTransferError,Plate_type,10)
        # Shouldn't be able to create a T-well plate, that makes no sense
        self.assertRaises(PlateTransferError,Plate_type,'T')


class Testing__generate_outfile_names(unittest.TestCase):
    """ Unit-tests for the generate_outfile_names function. """

    def test__single_Biomek_file(self):
        assert generate_outfile_names('X',0,0) == ('X.txt',['X_Biomek.csv'],[])

    def test__different_input_formats(self):
        # file_plate_names can be specified in multiple ways
        assert generate_outfile_names('X',1,0,1,['A']) == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,1,'A') == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],[])

    def test__mirroring(self):
        # mirroring
        assert generate_outfile_names('X',0,1) == ('X.txt',['X_Biomek.csv'],['X_Biomek_mirror.csv'])
        assert generate_outfile_names('X',1,1,1,['A']) == ('X.txt',['X_Biomek_A.csv'],['X_Biomek_mirror_A.csv'])
        assert generate_outfile_names('X',1,1,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                               ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                           ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],
                                                         ['X_Biomek_mirror_A1.csv','X_Biomek_mirror_A2.csv'])

    def test__second_arg_true_requires_last_two_args(self):
        # the last two args must be given if second arg is True
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0,['A'])
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1,['A'])


class Testing__get_plate_name_list_from_input(unittest.TestCase):
    """ Unit-tests for the get_plate_name_list_from_input function. """

    def test__input_types(self):
        # input can be a list of appropriate length (return unchanged) or a string (split on ,)
        assert get_plate_name_list_from_input(1,['A']) == ['A']
        assert get_plate_name_list_from_input(1,'A') == ['A']
        assert get_plate_name_list_from_input(4,['A','B','C','D']) == ['A','B','C','D']
        assert get_plate_name_list_from_input(4,'A,B,C,D') == ['A','B','C','D']
        # passing a one-element string is allowed for any N - plates are numbered automatically
        assert get_plate_name_list_from_input(4,'A') == ['A1','A2','A3','A4']
        # if the second arg is a LIST, not string, with one element, fail 
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,['A'])

    def test__correct_number_of_plate_names(self):
        # first arg should match length of second arg (unless the latter is a one-element string)
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 1,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 1,['A','B','C'])
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 2,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 2,['A','B','C'])
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,['A','B','C'])

    def test__duplicates_not_allowed(self):
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 3,'A,A,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 3,['A','A','C'])


class Testing__numbers_to_plate_and_well_IDs(unittest.TestCase):
    """ Unit-tests for the numbers_to_plate_and_well_IDs function. """

    def test__correct_cases(self):
        assert numbers_to_plate_and_well_IDs(10, 6, 2, ['plate1','plate2']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate2,A1', 'plate2,A2', 'plate2,A3', 'plate2,B1']
        assert numbers_to_plate_and_well_IDs(10, 24, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate1,B4']
        assert numbers_to_plate_and_well_IDs(10, 96, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,A7', 'plate1,A8', 'plate1,A9', 'plate1,A10']

    def test__bad_input(self):
        # Shouldn't work, N_plates doesn't match the plate ID list
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 10, 6, 2, ['plate1'])
        # Shouldn't work, not enough plates
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 20, 6, 2, ['plate1','plate2'])
        # Shouldn't work, too many plates
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 2, 6, 2, ['plate1','plate2'])
        # Shouldn't work, 10 isn't a valid plate size
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 2, 10, 2, ['plate1','plate2'])


class Testing__assign_codewords(unittest.TestCase):
    """ Unit-tests for the assign_codewords function. """

    def setUp(self):
        # make some test binary codes
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        [self.b01,self.b10,self.b11,self.b00] = [b01,b10,b11,b00]
        self.b_list_with_00 = [b01,b10,b11,b00]
        self.b_list_without_00 = [b01,b10,b11]
        self.B_with_00 = binary_code_utilities.Binary_code(2,self.b_list_with_00)
        self.B_without_00 = binary_code_utilities.Binary_code(2,self.b_list_without_00)
        [b01111,b10000,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['01111','10000','10001']]
        [self.b01111,self.b10000,self.b10001] = [b01111,b10000,b10001]
        self.b_list_longer = [b01111,b10000,b10001]
        self.B_longer = binary_code_utilities.Binary_code(5,self.b_list_longer)

    def test__result_is_subset_of_right_length(self):
        """ Result should be a subset of the expected set and of the expected length (don't check the ordering yet). """
        # (should work the same regardless of remove_low value)
        for r in [True,False]:
            s = set(assign_codewords(3,2,self.B_without_00,remove_low=r,quiet=True))
            assert s == set(self.b_list_without_00) and len(s) == 3
            s = set(assign_codewords(2,2,self.B_without_00,remove_low=r,quiet=True))
            assert s.issubset(set(self.b_list_without_00)) and len(s) == 2
            s = set(assign_codewords(1,2,self.B_without_00,remove_low=r,quiet=True))
            assert s.issubset(set(self.b_list_without_00)) and len(s) == 1
            s = set(assign_codewords(3,5,self.B_longer,remove_low=r,quiet=True))
            assert s == set(self.b_list_longer) and len(s) == 3
            s = set(assign_codewords(2,5,self.B_longer,remove_low=r,quiet=True))
            assert s.issubset(set(self.b_list_longer)) and len(s) == 2
            s = set(assign_codewords(1,5,self.B_longer,remove_low=r,quiet=True))
            assert s.issubset(set(self.b_list_longer)) and len(s) == 1

    def test__all_zero_codeword_always_removed(self):
        """ The all-zero codeword should always be thrown away. """
        # (should work the same regardless of remove_low value)
        for r in [True,False]:
            assert self.b00 not in set(assign_codewords(3,2,self.B_with_00,remove_low=r,quiet=True))
            assert self.b00 not in set(assign_codewords(2,2,self.B_with_00,remove_low=r,quiet=True))
            assert self.b00 not in set(assign_codewords(1,2,self.B_with_00,remove_low=r,quiet=True))

    def test__words_sorted_by_weight(self):
        """ The words should be sorted by weight, with the low-weight ones taken first (reverse if remove_low is True).
        (Assume the order of same-weight words can be arbitrary - accept both cases.) """
        assert assign_codewords(1,2,self.B_without_00,quiet=True) in [ [self.b01], [self.b10] ]
        assert assign_codewords(2,2,self.B_without_00,quiet=True) in [ [self.b01,self.b10], [self.b10,self.b01] ]
        assert assign_codewords(3,2,self.B_without_00,quiet=True) in [[self.b01,self.b10,self.b11], 
                                                                      [self.b10,self.b01,self.b11]]
        assert assign_codewords(1,5,self.B_longer,quiet=True) == [self.b10000]
        assert assign_codewords(2,5,self.B_longer,quiet=True) == [self.b10000,self.b10001]
        assert assign_codewords(3,5,self.B_longer,quiet=True) == [self.b10000,self.b10001,self.b01111]
        # if remove_low is set to True, the sorting by weight should be reversed.
        assert assign_codewords(1,2,self.B_without_00,remove_low=True,quiet=True) == [self.b11]
        assert assign_codewords(2,2,self.B_without_00,remove_low=True,quiet=True) in [[self.b11,self.b10], 
                                                                                      [self.b11,self.b01]]
        assert assign_codewords(3,2,self.B_without_00,remove_low=True,quiet=True) in [[self.b11,self.b01,self.b10], 
                                                                                      [self.b11,self.b10,self.b01]]
        assert assign_codewords(1,5,self.B_longer,remove_low=True,quiet=True) == [self.b01111]
        assert assign_codewords(2,5,self.B_longer,remove_low=True,quiet=True) == [self.b01111,self.b10001]
        assert assign_codewords(3,5,self.B_longer,remove_low=True,quiet=True) == [self.b01111,self.b10001,self.b10000]

    def test__same_weight_words_sorted_lexicographically(self):
        """ Words of the same weight should be sorted lexicographically (but I'm not sure we want to rely on that). """
        assert assign_codewords(3,2,self.B_without_00,quiet=True) == [self.b01,self.b10,self.b11]
        assert assign_codewords(3,2,self.B_without_00,remove_low=True,quiet=True) == [self.b11,self.b01,self.b10]

    def test__same_weight_words_sorted_lexicographically(self):
        """ Function should fail when the number of samples or pools doesn't match. """
        # (should work the same regardless of remove_low value)
        for r in [True,False]:
            # Should be too many samples for given code
            self.assertRaises(PlateTransferError, assign_codewords, 4,2,self.B_without_00,remove_low=r,quiet=True)
            # Should be too many samples after '00' removal
            self.assertRaises(PlateTransferError, assign_codewords, 4,2,self.B_with_00,remove_low=r,quiet=True)
            # Number of pools doesn't match code length
            self.assertRaises(PlateTransferError, assign_codewords, 4,1,self.B_without_00,remove_low=r,quiet=True)
            self.assertRaises(PlateTransferError, assign_codewords, 4,3,self.B_without_00,remove_low=r,quiet=True)


class Testing__make_Biomek_file_commands(unittest.TestCase):
    """ Unit-tests for the make_Biomek_file_commands function. """

    def test__basic_functionality(self):
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        # basic functionality for combinatorial pooling (note that 'x', 'A' etc here would really be 'plate1,A1' or such)
        assert make_Biomek_file_commands([b10],['x'],['A','B'],5) == ['x,A,5']
        assert make_Biomek_file_commands([b01],['x'],['A','B'],5) == ['x,B,5']
        assert make_Biomek_file_commands([b11],['x'],['A','B'],5) == ['x,A,5','x,B,5']
        assert make_Biomek_file_commands([b01,b10,b11],['x','y','z'],['A','B'],5) == ['x,B,5','y,A,5','z,A,5','z,B,5']
        assert make_Biomek_file_commands([b11,b10,b01],['x','y','z'],['A','B'],5) == ['x,A,5','x,B,5','y,A,5','z,B,5']
        assert make_Biomek_file_commands([b01,b01,b01],['x','y','z'],['A','B'],5) == ['x,B,5','y,B,5','z,B,5']

    def test__fail_for_length_mismatches(self):
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        [b1,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['1','10001']]
        # Lengths of codewords and sample_positions are mismatched
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b11],['x','y','z'],['A','B'],5)
        # Lengths of codewords and sample_positions are mismatched
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b11],['x','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b1],['x','y','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b10001],['x','y','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b11],['x','y','z'],['B'],5)


class Testing__split_command_list_by_source(unittest.TestCase):
    """ Unit-tests for the split_command_list_by_source function. """

    def test__empty_list(self):
        """ empty list -> empty dictionary """
        assert split_command_list_by_source([]) == {}

    def test__single_plate(self):
        """ if there's only one plate, the result should be a one-item dictionary. """
        assert split_command_list_by_source(['p1,A1,x,5','p1,A2,y,5']) == {'p1':['p1,A1,x,5','p1,A2,y,5']}

    def test__multiple_plates(self):
        """ if there are multiple plates, return one dict per plate. """
        assert split_command_list_by_source(['p1,A1,x,5','p2,A1,y,5']) == {'p1':['p1,A1,x,5'], 'p2': ['p2,A1,y,5']}
        assert split_command_list_by_source(['p1,A1,x,5','p2,A1,y,5','p2,A2,y,5']) == {'p1':['p1,A1,x,5'], 
                                                                                       'p2': ['p2,A1,y,5','p2,A2,y,5']}


class Testing__split_command_list_to_max_commands(unittest.TestCase):

    def test__bad_inputs(self):
        """ N must be a positive integer (yes, I'm not checking for other wrong types, sue me). """
        for N in [-1,0]:
            self.assertRaises(PlateTransferError, split_command_list_to_max_commands, [], N)
            self.assertRaises(PlateTransferError, split_command_list_to_max_commands, ['a','b','c'], N)

    def test__max_greater_than_len(self):
        """ if N<len(list), return the original list. """
        for N in range(1,10):
            assert split_command_list_to_max_commands([], N) == []
            assert split_command_list_to_max_commands(['a','b','c'], N+3) == [['a','b','c']]

    def test__various(self):
        """ Testing specific cases. """
        assert split_command_list_to_max_commands(['a','b','c'], 1) == [['a'],['b'],['c']]
        assert split_command_list_to_max_commands(['a','b','c'], 2) == [['a','b'],['c']]
        assert split_command_list_to_max_commands(['a','b','c'], 3) == [['a','b','c']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 1) == [['a'],['b'],['c'],['d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 2) == [['a','b'],['c','d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 3) == [['a','b'],['c','d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 4) == [['a','b','c','d']]

def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    if not os.access("./error-correcting_codes",os.F_OK):
        return "Error: there is not error-correcting_codes folder in this directory - can't run tests."
    if os.access("./test_outputs",os.F_OK):
        print("Test output files will be saved in the test_outputs directory (already present - removing it now).")
        shutil.rmtree("./test_outputs")
    else:
        print("Test output files will be saved in the test_outputs directory (not present - creating it now).")
    os.mkdir("./test_outputs")
    # MAYBE-TODO Allow the test_outputs folder name to be given as an argument/option?  How?

    test_runs = ["-n 63 -N 15 -P 3 -i Source1 -o -C error-correcting_codes/15-6-6_generator test_outputs/test1 -q", 
                 "-n 63 -N 15 -P 3 -i Source1 -o -M -C error-correcting_codes/15-6-6_generator test_outputs/test2 -q", 
                "-n 384 -N 18 -p 4 -P 3 -i Source -m -C error-correcting_codes/18-9-6_generator test_outputs/test3 -q"]
    # MAYBE-TODO add name/description strings to the test cases?
    parser = define_option_parser()
    for test_run in test_runs:
        print(" * New test run, with arguments: %s"%test_run)
        # regenerate options with test argument string
        (options, args) = parser.parse_args(test_run.split())
        run_main_function(parser,options,args)
    print("*** Test runs finished. If you didn't get any errors, that's good (warnings are all right). "
          + "Check the output files to make sure they look reasonable (this is NOT done automatically!).")
    return 0
    # TODO move this to a separate function!!
    # TODO write full tests that compare the outputs to expected files made by hand!!!
    # TODO I'd like to be able to compare files even if they have different date headers - I could use the command-line diff utility and its --ignore-matching-lines=<REGEX> option, or python difflib package (ndiff or SequenceMatcher or something http://docs.python.org/library/difflib.html) and the linejunk function.  Or just write something to remove the "date" lines from the output files before doing the comparison? No, that's bad.
    # MAYBE-TODO Set up a more complicated file-to-reference comparison myself - like compare each reference/output line normally (i.e. they have to be identical), UNLESS the reference file line starts with <REGEX>, in which case check if the output file line matches the regex given in the reference file line.  See my stackoverflow question http://stackoverflow.com/questions/9726214/testing-full-program-by-comparing-output-file-to-reference-file-whats-it-calle

### General functions (not input/output or optparse-related or testing or main), in approximate order of use

def generate_outfile_names(outfile_basename,if_multiple_files,if_mirror_files,number_of_files=None,file_plate_names=None):
    """ Given the base outfile name, generate full outfile names: (general_outfile, Biomek, Biomek_mirror).
    If if_multiple_files is false:  X -> (X.txt, [X_Biomek.csv], M)  (the last two arguments are ignored in this case)
    Otherwise, something like this:  X -> (X.txt, [X_Biomek_A.csv,X_Biomek_B.csv,...], M)  
      (assuming get_plate_name_list_from_input(number_of_files,file_plate_names) returns something like [A,B,...]) 
    If if_mirror_files is true, Biomek_mirror (M above) is []; otherwise it's the same as Biomek with a _mirror suffix."""
    ### 1) generate the basic outfile names
    outfile_data = outfile_basename+'.txt'
    ### 2) generate the basic outfile names
    if not if_multiple_files:
        # note: just ignore the last two arguments, they may be inconsistent with a single file, that's FINE. 
        outfiles_Biomek = [outfile_basename+'_Biomek.csv']
    else:
        if not (number_of_files and file_plate_names):
            raise PlateTransferError("If outputting multiple Biomek files, must specify number and names!")
        # file_plate_names here can still be a single string or any number of other things - get a list
        file_plate_names = get_plate_name_list_from_input(number_of_files, file_plate_names)
        outfiles_Biomek = [outfile_basename+'_Biomek_'+plate+'.csv' for plate in file_plate_names]
    ### 3) generate the mirror outfile names (if requested)
    if not if_mirror_files:
        outfiles_Biomek_mirror = []
    else:
        outfiles_Biomek_mirror = [name.replace('_Biomek','_Biomek_mirror') for name in outfiles_Biomek]
    return (outfile_data,outfiles_Biomek,outfiles_Biomek_mirror)


def get_plate_name_list_from_input(N_plates,ID_input):
    """ Return a list of plate names of length N_plates, generated using ID_input. 
    If ID_input is already a correct list, return it.  Otherwise assume it's a string: if splitting it 
    on commas yields N results, return those; if it yields a single result X, return [X1,X2,...,XN].
    Note that if N_plates is 1, a count won't be appended: (1,'X') returns ['X'], not ['X1']. """
    # if we already have a list with the right number of arguments, just return it after checking it makes sense
    if isinstance(ID_input,list): 
        if not len(ID_input)==len(set(ID_input)):
            raise PlateTransferError("You specified a list of plate names with duplicate values! Don't do that.")
        if not len(ID_input)==N_plates:
            raise PlateTransferError("Passed a plate name list of wrong length to get_plate_name_list_from_input!")
        return ID_input
    # otherwise generate list from string:
    split_IDs = ID_input.split(',')
    # either it's a comma-separated 'list' of plate names, of correct length, and should be treated as above
    #  (testing for len==N_plates first: (1,'X') returns ['X'], not ['X1'] - easier, since (1,['X']) returns ['X'] too.)
    if len(split_IDs) == N_plates:  
        if not len(split_IDs)==len(set(split_IDs)):
            raise PlateTransferError("You specified a list of plate names with duplicate values! Don't do that.")
        return split_IDs
    # or it's a single plate name and sequential numbers should be appended to it to get the list.
    elif len(split_IDs) == 1:
        return ['%s%d'%(ID_input,n+1) for n in range(N_plates)]
    else:
        raise PlateTransferError("Can't figure out how to name %s plates using input \"%s\"!"%(N_plates,ID_input))


def numbers_to_plate_and_well_IDs(N_samples, plate_size, N_plates, plate_IDs):
    """ Given the number of samples and plate information, return a list of with plate/well positions for each sample.
    The list will be of length N_samples, about like this: ['plate1,A1','plate1,A2',...,'plate1,H12','plate2,A1',...]. """

    if not len(plate_IDs)==N_plates:
        raise PlateTransferError("The number of plates must match the number of plate IDs provided!")
    if not plate_size in Plate_type.plate_sizes:
        raise PlateTransferError("Plate size must be an integer in %s!"%Plate_type.plate_sizes)
    if N_samples > plate_size*N_plates:
        raise PlateTransferError("Can't fit %s samples in %s %s-well plates!"%(N_samples,N_plates,plate_size))
    if N_samples <= plate_size*(N_plates-1):
        raise PlateTransferError("Why use %s %s-well plates "%(N_plates,plate_size)
                                     + "when you can fit %s samples in %s plates?"%(N_samples,N_plates-1))

    plate_type = plate_types[plate_size]
    position_list = []
    for i in range(N_samples):
        well_number = i % plate_size
        well_ID = plate_type.get_well_ID_from_number(well_number)
        plate_number = int(i/plate_size)
        plate_ID = plate_IDs[plate_number]
        position_list.append("%s,%s"%(plate_ID,well_ID))

    return position_list


def assign_codewords(N_samples, N_pools, binary_code, remove_low=False, quiet=False):
    """ Return a list of Binary_codeword objects, of length N_samples - one codeword per sample, ordered. """

    if not N_pools == binary_code.length:
        raise PlateTransferError("The codeword length in the provided binary code must be equal to N_pools!")
    if N_samples > binary_code.size():
        raise PlateTransferError("N_samples cannot exceed the size of the provided binary code!")
    if binary_code.remove_extreme_codeword(bit=0):
        if not quiet:
            print("Warning: the binary code provided contained the all-zero codeword, which cannot be used, "
                  + "as it would result in a sample being added to none of the pools and thus not sequenced. "
                  + "The all-zero codeword was removed.")
        if N_samples > binary_code.size():
            raise PlateTransferError("After removing the all-zero codeword, there aren't enough codewords for all "
                                     + "samples! Aborting. You could consider inverting the code to get around this.")
    if N_samples < binary_code.size() and not quiet:
        print("Warning: N_samples is lower than the size of the provided binary code - an arbitrary subset of codewords "
              + "will be used.  You may want to reduce your code size manually for improved attributes.")

    # get the desired number of samples from the binary code (returns a sorted list, original code is unchanged)
    #  (removing either the low-weight or high-weight depending on the remove_low argument value)
    codeword_list = binary_code.give_N_codeword_list_by_bit_sum(N_samples,remove_low=remove_low)
    return codeword_list


def make_Biomek_file_commands(sample_codewords, sample_positions, pool_positions, volume):
    """ Return a list of Biomek transfer commands to perform combinatorial pooling based on sample_codewords.

    Inputs:
     - sample_codewords - sequential list of Binary_codeword object corresponding to each sample
     - sample_positions - list of plate/well position strings (like "Source1,A4") for each sample, in the same order
     - pool_positions - same-format list of pool position strings
     - volume - integer giving the volume of all the transfers 

    Combinatorial pooling: the pools correspond to each bit of the codeword.  Sample A should be added to pool X
    whenever bit X of the codeword for sample A is 1. 
    Biomek command list format:  a list of strings of the form "plateA,wellA,plateX,wellX,volume" 
    where "plateA,wellA" is the value of sample_positions[A], and "plateX,wellX" is pool_positions[X]. """

    # make sure the inputs make sense
    if not len(sample_codewords)==len(sample_positions):
        raise PlateTransferError("The number of sample positions doesn't match the number of codewords!")
    if not set([len(x) for x in sample_codewords]) == set([len(pool_positions)]):
        raise PlateTransferError("Not all codeword lentgths match the number of pools *%s)!"%len(pool_positions))
        # note that the second set is always of size 1, so this implicitly makes sure all codewords are the same length

    Biomek_file_commands = []
    for (sample_number, (sample_codeword, sample_position)) in enumerate(zip(sample_codewords,sample_positions)):
        pools_to_add_sample_to = [pool_number for (pool_number,if_add) in enumerate(sample_codeword.list()) if if_add==1]
        for pool_number in pools_to_add_sample_to:
            Biomek_file_commands.append("%s,%s,%s"%(sample_position,pool_positions[pool_number],volume))
    # MAYBE-TODO might be nice to have info printed about the highest number of pools per sample and samples per pool, to help with volume calculations... And also the total number of transfers (maybe divided by 96) to figure out how many tip boxes will be needed
    return Biomek_file_commands
# TODO make this calculate and output two pieces of information: 
# - how many samples will be taken from each source position!
# - how many samples will end up in each destination position!
# Make the formatting useful-looking - possibly just give the range first.


def split_command_list_by_source(Biomek_file_commands):
    """ Split list of "x,_,_,_,_" strings into multiple lists by x value, return a (x_val: line_list) dictionary."""
    data_dict = defaultdict(lambda: [])
    for line in Biomek_file_commands:
        source_plate = line.split(',')[0]
        data_dict[source_plate].append(line)
    return data_dict


def split_command_list_to_max_commands(Biomek_file_commands, max_lines=300):
    """ Split list of strings into multiple lists no longer than max_lines; return list of those lists."""
    if max_lines<=0: raise PlateTransferError("max_lines must be a positive integer!")
    if len(Biomek_file_commands)==0:    return []
    N_lines = len(Biomek_file_commands)
    N_lists = int(ceil(float(N_lines)/max_lines))
    # rather than just take max_lines until nothing is left, split the original list up evenly: 
    #  even if max_lines is 3, len4 should become [len2,len2] rather than [len3,len1]
    N_lines_per_list = int(ceil(float(N_lines)/N_lists))
    #print "%s lines, %s max lines per list -> %s lists with <=%s lines"%(N_lines, max_lines, N_lists, N_lines_per_list)
    new_lists = []
    for i in range(N_lists):
        new_lists.append(Biomek_file_commands[i*N_lines_per_list : (i+1)*N_lines_per_list])
    return new_lists


### Input/output functions - no need/ability to unit-test, all the complicated functionality should be elsewhere.

def get_binary_code(length,listfile=None,matrixfile=None):
    """ Given a listfile or matrixfile name (but not both!), return the generated binary code."""
    if not (listfile or matrixfile):
        raise PlateTransferError("You must provide either a listfile or a matrixfile to generate a binary code!")
    if listfile and matrixfile:
        raise PlateTransferError("You cannot provide both a listfile and a matrixfile to generate a single binary code!")
    if listfile:
        method = 'listfile'
        infile = listfile
    elif matrixfile:
        method = 'matrixfile'
        infile = matrixfile
    binary_code = binary_code_utilities.Binary_code(length=length, val=infile, method=method)
    return binary_code


def write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, max_commands_per_file=0, Biomek_header=""):
    """ Write Biomek_file_commands to one or more Biomek outfiles with given header.  

    Biomek_file_commands and outfiles_Biomek must be lists of matching length (giving file contents and filenames).

    The outfiles_Biomek argument must be a list: containing a single element if there will be one outfile 
    (in which case Biomek_file_commands should be a single list of command strings), 
    or more for multiple outfiles (in which case Biomek_file_commands should be a dict of matching length, 
    with each key being the plate name which should match the Biomek file name, and each value being 
    a list of command strings to be written to the corresponding Biomek file). """
    # TODO finish rewriting docstring! Make sure it makes sense; add the max_commands_per_file bit.
    data_filename_list = []
    if len(outfiles_Biomek)==1:
        data_filename_list.append((Biomek_file_commands, outfiles_Biomek[0]))
    else:
        # split the commands by source
        transfer_file_command_sets = split_command_list_by_source(Biomek_file_commands)
        # Sort both the command set list (make list from dict first) and the Biomek outfile set.  
        #  (They're based on the same underlying plate names passed to the function, plus invariant prefixes/suffixes, 
        #  so they should always match once they're sorted, even if the plate names weren't sorted sensibly themselves.)
        transfer_file_command_sets = sorted(list(transfer_file_command_sets.items()))
        outfiles_Biomek.sort()
        # make sure the resulting lists match by length and by plate name
        if not len(transfer_file_command_sets) == len(outfiles_Biomek):
            raise PlateTransferError("The number of Biomek outfiles and command data sets provided doesn't match!")
        for ((set_name,_),outfilename) in zip(transfer_file_command_sets,outfiles_Biomek):
            if not set_name in outfilename:
                raise PlateTransferError("Can't match the Biomek file names and command sets!")
        # write each data set to the corresponding file
        for ((_,data),outfilename) in zip(transfer_file_command_sets,outfiles_Biomek):
            data_filename_list.append((data, outfilename))
    # now for each (data,filename) tuple, actually write it to a file (or to multiple files if splitting is necessary)
    for line_list,filename in data_filename_list:
        if max_commands_per_file==0:
            save_line_list_as_file(line_list, filename, header=Biomek_header)
        else:
            pass
            # TODO implement! Using split_command_list_to_max_commands; give a/b/c/ suffixes to outfiles.


def write_data_to_outfile(sample_codewords, sample_positions, pool_positions, 
                          outfile_data, outfiles_Biomek, options=None):
    """ Write data to main outfile: header information (command, path, date/time, options - all as #-start comments), 
    sample numbers/positions/codewords, and pool numbers/positions. """
    # print all the usual header information (command, path, date/time, options)
    OUTFILE = open(outfile_data,'w')
    write_header_data(OUTFILE,options)
    OUTFILE.write("# Corresponding Biomek file(s): %s\n"%outfiles_Biomek)
    # write the sample_number_position_codeword_list and pool_number_position_list data
    OUTFILE.write("\nsample_number,\tplate_and_well_position,\tcodeword\n")
    for (number,(codeword,position)) in enumerate(zip(sample_codewords,sample_positions)):
        OUTFILE.write("%s,\t%s,\t%s\n"%(number,position,codeword.string()))
    OUTFILE.write("\npool_number,\tplate_and_well_position\n")
    for (number,position) in enumerate(pool_positions):
        OUTFILE.write("%s,\t%s\n"%(number,position))
    # MAYBE-TODO print the mirror destination plate positions as well?
    OUTFILE.close()


### Option parser functions - no real need to unit-test, and it'd be complicated and weird.

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    parser.add_option('-q','--quiet', action='store_true', default=False, help="Don't print warnings (default False).")
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default False).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run with a set of predetermined realistic options to make sure it more or less works "
                      + "(default False). Ignores all other options/arguments; output files will start with 'test'.")

    parser.add_option('-n','--number_of_samples', type='int', metavar='N', help="Number of samples to pool (required).")
    parser.add_option('-N','--number_of_pools', type='int', metavar='M', help="Number of resulting pools (required).")
    parser.add_option('-c','--binary_code_list_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a list of codewords).")
    parser.add_option('-C','--binary_code_generator_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a generator matrix).")
    parser.add_option('-M','--add_mirror_pooling_files', action='store_true', default=False, 
                      help="In addition to the normal Biomek file, also make files with commands for a 'mirrored' set: " 
                      + "if sample A is in pool B in the normal set it isn't in the mirrored set, and vice versa.")
    parser.add_option('-u','--mirror_pool_plate_suffix', default='_mirror', metavar='S', 
                      help="Append S to the pool plate names for the mirror pooling Biomek files (default %default).")

    parser.add_option('-s','--size_of_sample_plates', type='int', default=96, metavar='M', 
                      help="Sample (source) plate size (from %s)"%Plate_type.plate_sizes + " (default %default)")
    parser.add_option('-S','--size_of_pool_plates', type='int', default=6, metavar='M', 
                      help="Pool (destination) plate size (from %s)"%Plate_type.plate_sizes + " (default %default)")
    parser.add_option('-p','--number_of_sample_plates', type='int', default=1, metavar='N', 
                      help="Total number of sample (source) plates to use (default %default).")
    parser.add_option('-P','--number_of_pool_plates', type='int', default=4, metavar='N', 
                      help="Total number of pool (destination) plates to use (default %default).")
    parser.add_option('-i','--sample_plate_IDs', default='Source', metavar='S', 
                      help="Sample plate IDs (must match the IDs in the Biomek deck setup): either a comma-separated list "
                      + "with the number of values matching -n ('A,B,C,D' - plates will be named A, B, C and D), "
                      + "or a single string with no commas ('x' - plates will be named x1, x2 etc.) (default %default).")
    parser.add_option('-I','--pool_plate_IDs', default='Destination', metavar='S', 
                      help="Pool plate IDs - see -i for details on how it works. Default %default.")

    parser.add_option('-m','--multiple_Biomek_files', action='store_true', default=True,
                      help="Generate multiple Biomek files, one per source plate (on by default).")
    parser.add_option('-o','--one_Biomek_file', action='store_false', dest='multiple_Biomek_files',
                      help="Generate a single Biomek file regardless of number of plates involved (off by default).")
    parser.add_option('-x','--max_commands_per_file', type='int', default=0, metavar='N',
                      help="Split Biomek command files so they contain to more than N commands "
                          +"(use a/b/c/... suffixes for the split files) (0 means no maximum; default %default).")

    parser.add_option('-v','--volume_per_transfer', type='int', default=20, metavar='V', 
                      help="Liquid volume to use for each sample-to-pool transfer (default %default).")
    parser.add_option('-H','--Biomek_file_header', default="SourcePlt,SourceWell,DestPlt,DestWell,Volume", metavar='S', 
                      help="Header line for the Biomek transfer file (won't be printed if set to ''; default %default).")

    # MAYBE-TODO add a minimum Hamming distance option, make the program check that all pairs in the final codeword set satisfy that?
    # MAYBE-TODO get histograms of Hamming distances and bit sums?  (at least as numbers, possibly plots?)

    # MAYBE-TODO implement more ways of dealing with the clonality issue?  
    #   1. Simple check - check that the bitwise OR of no two codewords generates a codeword that's in the set
    #   2. More complicated check - minimum Hamming distance we require between any codeword and any bitwise OR of two codewords (would have to be specified as an option).  Or just get a histogram of that printed. 
    #   3. Actually try to get a set that doesn't have clonality issues? (based on some minimum distance given as an option).  (The actual mechanism for this would have to be implemented in binary_code_utilities.py - reduce_by_Hamming_distance or something else, see possible methods in the notes above)

    return parser


def check_options_and_args(parser,options,args):
    """ Take optparse parser/options/args, check number of args, check required options and option conflicts. """
    # TODO rewrite this so it doesn't require the parser as an argument, by moving the first try/except to __main__?  Why is it a problem, really? 

    try:
        [outfile_basename] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: There must be exactly one output file base name (shown as X in the examples). "
                 + "The general output file will be X.txt, the Biomek output file will be X_Biomek.csv "
                 + "(or X_Biomek_Source1.csv, X_Biomek_Source2.csv, etc, with the -m option).")

    if not all([ (s in Plate_type.plate_sizes) for s in [options.size_of_sample_plates, options.size_of_pool_plates] ]):
        sys.exit("Plate sizes (-s and -S) must be in integer in %s! No other sizes are defined."%Plate_type.plate_sizes)
    if not (options.number_of_samples>0 and options.number_of_pools>0):
        sys.exit("Positive -n and -N values required!")
    if not bool(options.binary_code_list_file) ^ bool(options.binary_code_generator_file):  # is xor
        sys.exit("Exactly one of -c and -C must be provided!")

    # MAYBE-TODO could allow -p/-P to be automatically calculated from -n/-N and -s/-S?

    return options,outfile_basename


### Main functionality - difficult to unit-test due to optparse options as an argument, can just be tested with full runs.

def run_main_function(parser,options,args):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    Takes an optparse parser object, and the options object and arg list generated by the parser."""
    # MAYBE-TODO may be more convenient for interactive use if this just took an input string, and generated/defined the parser itself...
    options,outfile_basename = check_options_and_args(parser,options,args)
    outfiles = generate_outfile_names(outfile_basename, options.multiple_Biomek_files, options.add_mirror_pooling_files, 
                                      options.number_of_sample_plates, options.sample_plate_IDs)
    (outfile_data, outfiles_Biomek, outfiles_Biomek_mirror) = outfiles
    print("Output files: %s"%(outfiles,))
    # assign codewords to samples
    binary_code = get_binary_code(options.number_of_pools, 
                                  options.binary_code_list_file, options.binary_code_generator_file)
    sample_codewords = assign_codewords(options.number_of_samples, options.number_of_pools, binary_code, 
                                        quiet=options.quiet)
    # generate plate names from strings if they weren't given as lists
    input_plate_names = get_plate_name_list_from_input(options.number_of_sample_plates, options.sample_plate_IDs)
    output_plate_names = get_plate_name_list_from_input(options.number_of_pool_plates, options.pool_plate_IDs)
    # generate the plate+well position strings for each input sample and each output pool
    sample_positions = numbers_to_plate_and_well_IDs(options.number_of_samples, options.size_of_sample_plates, 
                                                     options.number_of_sample_plates, input_plate_names)
    pool_positions = numbers_to_plate_and_well_IDs(options.number_of_pools, options.size_of_pool_plates, 
                                                   options.number_of_pool_plates, output_plate_names)
    # generate the Biomek transfer command list based on sample codewords and sample/pool positions
    Biomek_file_commands = make_Biomek_file_commands(sample_codewords, sample_positions, 
                                                                            pool_positions, options.volume_per_transfer)
    # write to outfiles
    write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, 
                               options.max_commands_per_file, options.Biomek_file_header)
    write_data_to_outfile(sample_codewords, sample_positions, pool_positions, 
                          outfile_data, outfiles_Biomek+outfiles_Biomek_mirror, options)
    # generate mirror Biomek files: invert the codewords, add suffix to pool plate names, run same functions.
    if options.add_mirror_pooling_files:
        mirror_sample_codewords = [~codeword for codeword in sample_codewords]
        mirror_output_plate_names = [plate_name+options.mirror_pool_plate_suffix for plate_name in output_plate_names]
        mirror_pool_positions = numbers_to_plate_and_well_IDs(options.number_of_pools, options.size_of_pool_plates, 
                                                              options.number_of_pool_plates, mirror_output_plate_names)
        mirror_Biomek_file_commands = make_Biomek_file_commands(mirror_sample_codewords, sample_positions, 
                                                                mirror_pool_positions, options.volume_per_transfer)
        write_data_to_Biomek_files(outfiles_Biomek_mirror, mirror_Biomek_file_commands, 
                                   options.max_commands_per_file, options.Biomek_file_header)


if __name__=='__main__':

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # Unit-testing - don't even look for more options/arguments, just run the test suite
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments (including -T), "
              + "running the built-in simple test suite.")
        print("Defined plate sizes: %s"%Plate_type.plate_sizes)
        # tun unittest.main, passing it no arguments (by default it takes sys.argv and complains about the -t)
        unittest.main(argv=[sys.argv[0]])
        # MAYBE-TODO unittest.main automatically quits - there's no way of doing -t and -T at once.  Do I care?
        #   may be fixed in a future version, and there is a patch: http://bugs.python.org/issue3379


    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs.")
        test_result = do_test_run()
        sys.exit(test_result)

    # If it's not a test run, just run the main functionality
    run_main_function(parser,options,args)

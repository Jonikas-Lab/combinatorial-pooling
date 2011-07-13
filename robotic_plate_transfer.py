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

import sys
from collections import defaultdict
import binary_code_utilities
from general_utilities import invert_list_to_dict, save_line_list_as_file, write_header_data

from string import ascii_uppercase as letters    # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

class PlateTransferError(Exception):
    """ Exception for this file (does nothing interesting)."""
    pass

plate_types_rows_columns = {6: (2,3), 24: (4,6), 96: (8,12), 384: (16,24)}
plate_sizes = sorted(plate_types_rows_columns.keys())


class Plate_type:
    """ A plate type (384-well, 96-well, 24-well or 6-well)."""

    def __init__(self,size):
        """ Go from the size to row and column numbers, call make_well_ID_list_and_dict."""
        self.size = size
        try:                
            self.rows, self.columns = plate_types_rows_columns[size]
        except KeyError:    
            raise PlateTransferError("Plate size must be in integer in %s!"%plate_sizes)
        assert self.rows*self.columns == self.size
        self.make_well_ID_list_and_dict()
    
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
plate_types = dict([(n,Plate_type(n)) for n in plate_sizes])


### Unit-tests for classes/setup above and "general functions" below

def test_functionality():
    import testing_utilities

    if True:
        print("Testing Plate_type class...")
        ### all the plate types have been generated already on module import, so just test them
        # try the first, second, first-in-second-row, and last well of each plate size
        if True:
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
        # and in reverse
        if True:
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
        # getting a too-high number should fail
        if True:
            testing_utilities.call_should_fail(plate_types[6].get_well_ID_from_number, [6], PlateTransferError, 
                                               message="Shouldn't be able to get well 6 (0-based) from 6-well plate!")
            testing_utilities.call_should_fail(plate_types[24].get_well_ID_from_number, [24], PlateTransferError, 
                                               message="Shouldn't be able to get well 24 (0-based) from 24-well plate!")
            testing_utilities.call_should_fail(plate_types[96].get_well_ID_from_number, [96], PlateTransferError, 
                                               message="Shouldn't be able to get well 96 (0-based) from 96-well plate!")
            testing_utilities.call_should_fail(plate_types[384].get_well_ID_from_number, [384], PlateTransferError, 
                                               message="Shouldn't be able to get well 384 (0-based) from 384-well plate!")

        # other things that should fail - making a new plate type with a weird size
        if True:
            testing_utilities.call_should_fail(Plate_type,[0],PlateTransferError,
                                               message="Shouldn't ever be able to create a 0-well plate!")
            testing_utilities.call_should_fail(Plate_type,[10],PlateTransferError,
                                               message="Shouldn't be able to create a 10-well plate, it wasn't defined!")
            testing_utilities.call_should_fail(Plate_type,['T'],PlateTransferError,
                                               message="Shouldn't be able to create a T-well plate, that makes no sense!")
        print("...DONE")

    if True:
        print("Testing generate_outfile_names function...")
        # single Biomek file - easy
        assert generate_outfile_names('X',0,0) == ('X.txt',['X_Biomek.csv'],[])
        # file_plate_names can be specified in multiple ways
        assert generate_outfile_names('X',1,0,1,['A']) == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,1,'A') == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],[])
        # mirroring
        assert generate_outfile_names('X',0,1) == ('X.txt',['X_Biomek.csv'],['X_Biomek_mirror.csv'])
        assert generate_outfile_names('X',1,1,1,['A']) == ('X.txt',['X_Biomek_A.csv'],['X_Biomek_mirror_A.csv'])
        assert generate_outfile_names('X',1,1,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                               ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                           ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],
                                                         ['X_Biomek_mirror_A1.csv','X_Biomek_mirror_A2.csv'])
        # the last two args must be given if second arg is True
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,0),PlateTransferError)
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,0,1),PlateTransferError)
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,0,['A']),PlateTransferError)
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,1),PlateTransferError)
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,1,1),PlateTransferError)
        testing_utilities.call_should_fail(generate_outfile_names,('X',1,1,['A']),PlateTransferError)
        print("...DONE")

    if True:
        print("Testing get_plate_name_list_from_input function...")
        # input can be a list of appropriate length (return unchanged) or a string (split on ,)
        assert get_plate_name_list_from_input(1,['A']) == ['A']
        assert get_plate_name_list_from_input(1,'A') == ['A']
        assert get_plate_name_list_from_input(4,['A','B','C','D']) == ['A','B','C','D']
        assert get_plate_name_list_from_input(4,'A,B,C,D') == ['A','B','C','D']
        assert get_plate_name_list_from_input(4,'A') == ['A1','A2','A3','A4']
        # first arg should match length of second arg (unless the latter is 1)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(1,'A,B,C'),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(1,['A','B','C']),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(2,'A,B,C'),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(2,['A','B','C']),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(4,'A,B,C'),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(4,['A','B','C']),PlateTransferError)
        # duplicate values not allowed
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(3,'A,A,C'),PlateTransferError)
        testing_utilities.call_should_fail(get_plate_name_list_from_input,(3,['A','A','C']),PlateTransferError)
        print("...DONE")

    if True:
        print("Testing numbers_to_plate_and_well_IDs function...")
        assert numbers_to_plate_and_well_IDs(10, 6, 2, ['plate1','plate2']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate2,A1', 'plate2,A2', 'plate2,A3', 'plate2,B1']
        assert numbers_to_plate_and_well_IDs(10, 24, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate1,B4']
        assert numbers_to_plate_and_well_IDs(10, 96, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,A7', 'plate1,A8', 'plate1,A9', 'plate1,A10']
        # things that should fail for various reasons
        testing_utilities.call_should_fail(numbers_to_plate_and_well_IDs, (10, 6, 2, ['plate1']), PlateTransferError, 
                                           message="Shouldn't work, N_plates doesn't match the plate ID list!")
        testing_utilities.call_should_fail(numbers_to_plate_and_well_IDs, (20, 6, 2, ['plate1','plate2']), 
                                           PlateTransferError, message="Shouldn't work, not enough plates!")
        testing_utilities.call_should_fail(numbers_to_plate_and_well_IDs, (2, 6, 2, ['plate1','plate2']), 
                                           PlateTransferError, message="Shouldn't work, too many plates!")
        testing_utilities.call_should_fail(numbers_to_plate_and_well_IDs, (2, 10, 2, ['plate1','plate2']), 
                                           PlateTransferError, message="Shouldn't work, 10 isn't a valid plate size!")
        print("...DONE")

    if True:
        print("Testing assign_codewords function... (ignore any warnings)")
        # make some text binary codes
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        b_list_with_00 = [b01,b10,b11,b00]
        b_list_without_00 = [b01,b10,b11]
        B_with_00 = binary_code_utilities.Binary_code(2,b_list_with_00)
        B_without_00 = binary_code_utilities.Binary_code(2,b_list_without_00)
        [b01111,b10000,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['01111','10000','10001']]
        b_list_longer = [b01111,b10000,b10001]
        B_longer = binary_code_utilities.Binary_code(5,b_list_longer)
        # first check that the result is a subset of the expected set - don't check the ordering
        assert set(assign_codewords(3,2,B_without_00)) == set(b_list_without_00)
        assert set(assign_codewords(2,2,B_without_00)).issubset(set(b_list_without_00))
        assert set(assign_codewords(1,2,B_without_00)).issubset(set(b_list_without_00))
        assert set(assign_codewords(3,5,B_longer)) == set(b_list_longer)
        assert set(assign_codewords(2,5,B_longer)).issubset(set(b_list_longer))
        assert set(assign_codewords(1,5,B_longer)).issubset(set(b_list_longer))
        # the all-zero codeword should always be thrown away
        assert set(assign_codewords(3,2,B_with_00)) == set(b_list_without_00)
        assert set(assign_codewords(2,2,B_with_00)).issubset(set(b_list_without_00))
        assert set(assign_codewords(1,2,B_with_00)).issubset(set(b_list_without_00))
        # the words are sorted first lexicographically, and then by weight (number of 1s)
        assert assign_codewords(3,2,B_without_00) == [b01,b10,b11]
        assert assign_codewords(2,2,B_without_00) == [b01,b10]
        assert assign_codewords(1,2,B_without_00) == [b01]
        assert assign_codewords(3,5,B_longer) == [b10000,b10001,b01111]
        assert assign_codewords(2,5,B_longer) == [b10000,b10001]
        assert assign_codewords(1,5,B_longer) == [b10000]
        # should fail when the number of samples or pools doesn't match
        testing_utilities.call_should_fail(assign_codewords,(4,2,B_without_00), PlateTransferError, 
                                           message="Should be too many samples for given code!")
        testing_utilities.call_should_fail(assign_codewords,(4,2,B_with_00), PlateTransferError, 
                                           message="Should be too many samples after '00' removal!")
        testing_utilities.call_should_fail(assign_codewords,(4,1,B_without_00), PlateTransferError, 
                                           message="Number of pools doesn't match code length!")
        testing_utilities.call_should_fail(assign_codewords,(4,3,B_without_00), PlateTransferError, 
                                           message="Number of pools doesn't match code length!")
        print("...DONE")

    if True:   
        print("Testing make_Biomek_file_commands function...")
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        # basic functionality for combinatorial pooling (note that 'x', 'A' etc here would really be 'plate1,A1' or such)
        assert make_Biomek_file_commands([b10],['x'],['A','B'],5) == ['x,A,5']
        assert make_Biomek_file_commands([b01],['x'],['A','B'],5) == ['x,B,5']
        assert make_Biomek_file_commands([b11],['x'],['A','B'],5) == ['x,A,5','x,B,5']
        assert make_Biomek_file_commands([b01,b10,b11],['x','y','z'],['A','B'],5) == ['x,B,5','y,A,5','z,A,5','z,B,5']
        assert make_Biomek_file_commands([b11,b10,b01],['x','y','z'],['A','B'],5) == ['x,A,5','x,B,5','y,A,5','z,B,5']
        assert make_Biomek_file_commands([b01,b01,b01],['x','y','z'],['A','B'],5) == ['x,B,5','y,B,5','z,B,5']
        # should fail if there are length mismatches between samples/codewords/pools
        [b1,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['1','10001']]
        testing_utilities.call_should_fail(make_Biomek_file_commands,([b01,b11],['x','y','z'],['A','B'],5), 
                                           PlateTransferError, message="Len of codewords and sample_positions mismatch!")
        testing_utilities.call_should_fail(make_Biomek_file_commands,([b01,b10,b11],['x','z'],['A','B'],5), 
                                           PlateTransferError, message="Len of codewords and sample_positions mismatch!")
        testing_utilities.call_should_fail(make_Biomek_file_commands,([b01,b10,b1],['x','y','z'],['A','B'],5), 
                                           PlateTransferError, message="Number of pools and codeword len mismatch!")
        testing_utilities.call_should_fail(make_Biomek_file_commands,([b01,b10,b10001],['x','y','z'],['A','B'],5), 
                                           PlateTransferError, message="Number of pools and codeword len mismatch!")
        testing_utilities.call_should_fail(make_Biomek_file_commands,([b01,b10,b11],['x','y','z'],['B'],5), 
                                           PlateTransferError, message="Number of pools and codeword len mismatch!")
        print("...DONE")

    if True:
        print("Testing split_command_list_by_source function...")
        # if there's only one plate, the result should be a one-item dictionary
        assert split_command_list_by_source(['p1,A1,x,5','p1,A2,y,5']) == {'p1':['p1,A1,x,5','p1,A2,y,5']}
        # if there are multiple plates, return one dict per plate
        assert split_command_list_by_source(['p1,A1,x,5','p2,A1,y,5']) == {'p1':['p1,A1,x,5'], 'p2': ['p2,A1,y,5']}
        print("...DONE")



### General functions (not input/output or optparse-related or testing or main), in approximate order of use

def generate_outfile_names(outfile_basename,if_multiple_files,if_mirror_files,number_of_files=None,file_plate_names=None):
    """ Given the base outfile name, generate full outfile names: (general_outfile, Biomek, Biomek_mirror).
    If if_multiple_files is false:  X -> (X.txt, [X_Biomek.csv], M)  (the last two arguments are ignored in this case)
    Otherwise, something like this:  X -> (X.txt, [X_Biomek_A.csv,X_Biomek_B.csv,...], M)  
      (assuming get_plate_name_list_from_input(number_of_files,file_plate_names) returns something like [A,B,...]) 
    If if_mirror_files is true, Biomek_mirror (M above) is []; otherwise it's the same as Biomek with a _mirror suffix."""

    outfile_data = outfile_basename+'.txt'
    if not if_multiple_files:
        # note: just ignore the last two arguments, they may be inconsistent with a single file, that's FINE. 
        outfiles_Biomek = [outfile_basename+'_Biomek.csv']
    else:
        if not (number_of_files and file_plate_names):
            raise PlateTransferError("If outputting multiple Biomek files, must specify number and names!")
        # file_plate_names here can still be a single string or any number of other things - get a list
        file_plate_names = get_plate_name_list_from_input(number_of_files, file_plate_names)
        outfiles_Biomek = [outfile_basename+'_Biomek_'+plate+'.csv' for plate in file_plate_names]
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
    if not plate_size in plate_sizes:
        raise PlateTransferError("Plate size must be an integer in %s!"%plate_sizes)
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


def assign_codewords(N_samples, N_pools, binary_code):
    """ Return a list of Binary_codeword objects, of length N_samples - one codeword per sample, ordered. """

    if not N_pools == binary_code.length:
        raise PlateTransferError("The codeword length in the provided binary code must be equal to N_pools!")
    if N_samples > binary_code.size():
        raise PlateTransferError("N_samples cannot exceed the size of the provided binary code!")
    if binary_code.remove_all_zero_codeword():
        print("Warning: the binary code provided contained the all-zero codeword, which cannot be used, "
              + "as it would result in a sample being added to none of the pools and thus not sequenced. "
              + "The all-zero codeword was removed.")
        if N_samples > binary_code.size():
            raise PlateTransferError("After removing the all-zero codeword, there aren't enough codewords for all "
                                     + "samples! Aborting. You could consider inverting the code to get around this.")
    if N_samples < binary_code.size():
        print("Warning: N_samples is lower than the size of the provided binary code - an arbitrary subset of codewords "
              + "will be used.  You may want to reduce your code size manually for improved attributes.")

    # make the codewords into a sorted list, to guarantee the result is always the same
    codewords = sorted(binary_code.codewords)
    # we may as well do something useful here - remove the codewords with highest weight, to avoid issues like the all-one codeword and generally to limit the number of samples per pool.  MAYBE-TODO pick another way?  Make it an option?
    codewords.sort(key = lambda codeword: codeword.weight())
    # pick the required number of codewords
    sample_codewords = codewords[:N_samples]
    return sample_codewords


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
    return Biomek_file_commands


def split_command_list_by_source(Biomek_file_commands):
    """ Split list of "x,_,_,_,_" strings into multiple lists by x value, return a (x_val: line_list) dictionary."""
    data_dict = defaultdict(lambda: [])
    for line in Biomek_file_commands:
        source_plate = line.split(',')[0]
        data_dict[source_plate].append(line)
    return data_dict


### Input/output functions - no need/ability to unit-test, all the complicated functionality should be elsewhere.

def get_binary_code(listfile=None,matrixfile=None):
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
    binary_code = binary_code_utilities.Binary_code(length=options.number_of_pools, val=infile, method=method)
    return binary_code


def write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, Biomek_header=""):
    """ Write Biomek_file_commands to one or more Biomek outfiles with given header.  
    The outfiles_Biomek argument must be a list: containing a single element if there will be one outfile 
    (in which case Biomek_file_commands should be a single list of command strings), 
    or more for multiple outfiles (in which case Biomek_file_commands should be a dict of matching length, 
    with each key being the plate name which should match the Biomek file name, and each value being 
    a list of command strings to be written to the corresponding Biomek file). """
    if len(outfiles_Biomek)==1:
        save_line_list_as_file(Biomek_file_commands, outfiles_Biomek[0], header=Biomek_header)
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
            save_line_list_as_file(data, outfilename, header=Biomek_header)


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
    OUTFILE.close()


### Option parser functions - no real need to unit-test, and it'd be complicated and weird.

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

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
    parser.add_option('-u','--mirror_pool_plate_suffix', default='', metavar='S', 
                      help="Append S to the pool plate names for the mirror pooling Biomek files (default empty).")

    parser.add_option('-s','--size_of_sample_plates', type='int', default=96, metavar='M', 
                      help="Sample (source) plate size (from %s)"%plate_sizes + " (default %default)")
    parser.add_option('-S','--size_of_pool_plates', type='int', default=6, metavar='M', 
                      help="Pool (destination) plate size (from %s)"%plate_sizes + " (default %default)")
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

    parser.add_option('-v','--volume_per_transfer', type='int', default=20, metavar='V', 
                      help="Liquid volume to use for each sample-to-pool transfer (default %default).")
    parser.add_option('-H','--Biomek_file_header', default="SourcePlt,SourceWell,DestPlt,DestWell,Volume", metavar='S', 
                      help="Header line for the Biomek transfer file (won't be printed if set to ''; default %default).")

    return parser


def check_options_and_args(parser,options,args):
    """ Take optparse parser/options/args, check number of args, check required options and option conflicts. """

    try:
        [outfile_basename] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: There must be exactly one output file base name (shown as X in the examples). "
                 + "The general output file will be X.txt, the Biomek output file will be X_Biomek.csv "
                 + "(or X_Biomek_Source1.csv, X_Biomek_Source2.csv, etc, with the -m option).")

    if options.size_of_sample_plates not in plate_sizes or options.size_of_pool_plates not in plate_sizes:
        parser.error("Plate sizes (-s and -S) must be in integer in %s! No other plate sizes are defined."%plate_sizes)
    if not (options.number_of_samples>0 and options.number_of_pools>0):
        parser.error("Positive -n and -N values required!")
    if not bool(options.binary_code_list_file) ^ bool(options.binary_code_generator_file):  # is xor
        parser.error("Exactly one of -c and -C must be provided!")
    if options.mirror_pool_plate_suffix and not options.add_mirror_pooling_files:
        parser.error("There's no point in specifying -u when -M isn't set!")

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
    sample_codewords = assign_codewords(options.number_of_samples, options.number_of_pools, 
                                       get_binary_code(options.binary_code_list_file, options.binary_code_generator_file))
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
    write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, options.Biomek_file_header)
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
        write_data_to_Biomek_files(outfiles_Biomek_mirror, mirror_Biomek_file_commands, options.Biomek_file_header)


if __name__=='__main__':

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # Unit-testing - don't even look for more options/arguments, just run the test suite
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in simple test suite.")
        print("Defined plate sizes: %s"%plate_sizes)
        test_functionality()
        # if not doing a test run, exit now; doing a normal run after unit-testing isn't allowed, but doing a test run is.
        if not options.test_run:
            sys.exit(0)

    # Test run: robotic_plate_transfer.py -n 63 -N 15 -P 3 -i Source1 -C error-correcting_codes/15-6-6_generator test
    # since we'll be redoing the parsing on an example string, the option value will be overwritten, so save it separately
    test_run = options.test_run     
    test_run_inputs = ["-n 63 -N 15 -P 3 -i Source1 -o -C error-correcting_codes/15-6-6_generator test1", 
                       "-n 63 -N 15 -P 3 -i Source1 -o -M -u _mirror -C error-correcting_codes/15-6-6_generator test2", 
                       "-n 384 -N 18 -p 4 -P 3 -i Source -m -C error-correcting_codes/18-9-6_generator test3"]
    # MAYBE-TODO the -C option values above will only work if we're in the directory where the script is - fix that?
    # MAYBE-TODO add name/description strings to the test cases?
    if test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs.")
        for test_input in test_run_inputs:
            print(" * New test run, with arguments: %s"%test_input)
            # regenerate options with test argument string
            (options, args) = parser.parse_args(test_input.split())
            run_main_function(parser,options,args)
        print("*** Test runs finished. If you didn't get any errors, that's good (warnings are all right). "
              + "Check the output files to make sure.")
        # MAYBE-TODO add a -q option to silence the warnings for testing?
        sys.exit(0)

    # If it's not a test run, just run the main functionality
    run_main_function(parser,options,args)

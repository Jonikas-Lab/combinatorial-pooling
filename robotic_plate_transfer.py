#!/usr/bin/env python2
"""

 -- Weronika Patena, Jonikas lab, Carnegie Institution, July 2011
"""
# TODO write documentation based on the main function!

import collections
import string
import sys

class PlateTransferError(Exception):
    """ Exception for this file (does nothing interesting)."""
    pass

letters = string.ascii_uppercase  # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

plate_types_rows_columns = {6: (2,3), 24: (4,6), 96: (8,12), 384: (16,24)}
plate_sizes = sorted(plate_types_rows_columns.keys())

print "Defined plate sizes: %s"%plate_sizes

def invert_list_to_dict(input_list):
    """ Given a list with no duplicates, return a dict mapping the values to list positions ([a,b,c] -> {a:1,b:2,c:3})."""
    if not len(set(input_list)) == len(input_list):
        raise PlateTransferError("Can't reliably invert a list with duplicate elements!")
    return dict([(value,index) for (index,value) in enumerate(input_list)])


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


# TODO now we want a function that takes a code (described how? the generator filename, or just [n,k,d] and the program will look for a file?), a number of input samples and a number of output pools, and the input/output plate types/counts/names, and a volume, and makes a Biomek transfer file!
# (TODO will also need to implement reduce_to_number on Binary_code, for when the code size is 1024 but we only want 700 or something - how should the codewords to use be picked?  By lower or higher weight, or whichever weight is more extreme?  Also remember to remove the all-zero codeword, or invert the code to get rid of it if we need exactly 2**k samples!)


def test():
    import testing_utilities

    # TODO should test invert_list_to_dict as well!
    
    print "Testing Plate_type class..."
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
    print "...DONE"

    if True:
        print "Testing numbers_to_plate_and_well_IDs function..."
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
        print "...DONE"

    # TODO add tests of new functions/classes when I add them!


if __name__=='__main__':
    # TODO once I finish implementing other functions here, make this actually do something useful when called from the command-line, and only do the testing if run with a single "test" argument or something!
    test()

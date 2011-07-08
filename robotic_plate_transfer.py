#!/usr/bin/env python2
"""
Given a number of samples, a number of pools, a binary code, and information about the physical plates used for the samples and pools, assign one codeword from the code to each sample.  Output two files:  a general output file containing all the script information and the plate/well positions and codewords for each sample (and plate/well positions for each pool), and a Biomek robot transfer file that when given the appropriate plates will execute the combinatorial pooling.

For more details see help for the samples_and_code_to_Biomek_file function in this file and the option help descriptions.

 -- Weronika Patena, Jonikas lab, Carnegie Institution, July 2011

USAGE:  robotic_plate_transfer.py [options] outfile_base_name
        robotic_plate_transfer.py -t
"""

import sys
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


def samples_and_code_to_Biomek_file(N_samples, N_pools, binary_code, volume,
                                    input_plate_N, input_plate_size, input_plate_names, 
                                    output_plate_N, output_plate_size, output_plate_names):
    """ Make a Biomek transfer file to combine N samples into M pools, using a provided binary code.
         
        Assumes the samples will be in sequential order over the provided number of same-sized input plates: 
        [plate1 well A1, plate1 well A2, ..., plate1 well B2, .., plate1 well A1, ...].
        The pools will be similarly distributed over the provided number of output plates.  The number and size of 
        output plates does not need to have any particular relationship to the number and size of input plates.

        The two _plate_names arguments can be provided as a list with the appropriate number of elements, 
        as a string with the appropriae number of elements separated by commas, or as a single-word string 
        (in which case sequential numbers will be appended to it to generate the plate names).

        Arbitrarily assigns codewords from the code to each sample; the pools correspond to each bit of the codeword. 
        Whenever a bit of the codeword for the sample is 1, that sample should be added to the corresponding pool, 
        so a line of the form "sample_plate_name,sample_plate_well,pool_plate_name,pool_plate_well,volume" is added to 
        the Biomek command list. 

        Return the Biomek transfer file command list, a list of (input_sample_number,plate_and_well_position_string,
        codeword) tuples to keep track of which sample corresponds to which codeword, and a list of (pool_number,
        plate_and_well_position_string) tuples to keep track of pool positions (the last is not strictly necessary, 
        since the pools will always be in sequential order over the plates/wells, but it can be helpful to have a list 
        instead of regenerating it every time). """

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

    # generate plate names from strings if they weren't given as lists
    if not len(input_plate_names)==input_plate_N:
        input_plate_names = get_plate_names_from_option_value(input_plate_N, input_plate_names)
    if not len(output_plate_names)==output_plate_N:
        output_plate_names = get_plate_names_from_option_value(output_plate_N, output_plate_names)

    # generate the plate+well position strings for each input sample and each output pool
    sample_positions = numbers_to_plate_and_well_IDs(N_samples, input_plate_size, input_plate_N, input_plate_names)
    pool_positions = numbers_to_plate_and_well_IDs(N_pools, output_plate_size, output_plate_N, output_plate_names)
    # keep track of sample - position - codeword correspondences!  Keep a list of tuples, print to file at the end.
    sample_number_position_codeword_list = zip(range(N_samples), sample_positions, sorted(list(binary_code.codewords)))
    transfer_file_command_list = []
    for (sample_number, sample_position, sample_codeword) in sample_number_position_codeword_list:
        pools_to_add_sample_to = [pool_number for (pool_number,if_add) in enumerate(sample_codeword.list()) if if_add==1]
        for pool_number in pools_to_add_sample_to:
            transfer_file_command_list.append("%s,%s,%s"%(sample_position,pool_positions[pool_number],volume))
    # make a list of output pools and positions to return
    pool_number_position_list = list(enumerate(pool_positions))
    # for the return, we want the codewords to be strings, not objects
    sample_number_position_codeword_list = [(n,p,c.string()) for (n,p,c) in sample_number_position_codeword_list]
    return transfer_file_command_list, sample_number_position_codeword_list, pool_number_position_list


def get_plate_names_from_option_value(N_plates,ID_string):
    """ Return a list of length N_plates generated using ID_string: if splitting ID_string on commas yields N results, 
    return those; if it yields a single result X, return [X1,X2,...,XN]."""
    split_IDs = ID_string.split(',')
    if len(split_IDs) == N_plates:  
        return split_IDs
    elif len(split_IDs) == 1:
        return ['%s%d'%(ID_string,n+1) for n in range(N_plates)]
    else:
        raise PlateTransferError("Can't figure out how to name %s plates using string \"%s\"!"%(N_plates,ID_string))


def test_functionality():
    import testing_utilities

    if True:
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

    if True:
        print "Testing samples_and_code_to_Biomek_file function... (ignore any warnings)"
        import binary_code_utilities
        B1 = binary_code_utilities.Binary_code(2,['01','10','11','00'])
        command_list, sample_list, pool_list = samples_and_code_to_Biomek_file(3,2,B1,20, 1,24,['Src1'], 1,6,['Dest1'])
        assert command_list == ['Src1,A1,Dest1,A2,20', 'Src1,A2,Dest1,A1,20', 'Src1,A3,Dest1,A1,20', 'Src1,A3,Dest1,A2,20']
        assert sample_list == [(0, 'Src1,A1', '01'), (1, 'Src1,A2', '10'), (2, 'Src1,A3', '11')]
        assert pool_list == [(0, 'Dest1,A1'), (1, 'Dest1,A2')]
        testing_utilities.call_should_fail(samples_and_code_to_Biomek_file,(5,2,B1,20, 1,24,['Src1'], 1,6,['Dest1']),
                                           PlateTransferError, message="Should be too many samples for given code!")
        testing_utilities.call_should_fail(samples_and_code_to_Biomek_file,(4,2,B1,20, 1,24,['Src1'], 1,6,['Dest1']),
                                           PlateTransferError, message="Should be too many samples after '00' removal!")
        testing_utilities.call_should_fail(samples_and_code_to_Biomek_file,(3,1,B1,20, 1,24,['Src1'], 1,6,['Dest1']),
                                           PlateTransferError, message="Number of pools doesn't match code length!")
        testing_utilities.call_should_fail(samples_and_code_to_Biomek_file,(3,3,B1,20, 1,24,['Src1'], 1,6,['Dest1']),
                                           PlateTransferError, message="Number of pools doesn't match code length!")

        print "...DONE"

    if True:
        print "Testing get_plate_names_from_option_value function..."
        assert get_plate_names_from_option_value(4,'A,B,C,D') == ['A','B','C','D']
        assert get_plate_names_from_option_value(4,'s') == ['s1','s2','s3','s4']
        testing_utilities.call_should_fail(get_plate_names_from_option_value,(4,'A,B,C'),PlateTransferError)
        print "...DONE"


def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default False).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run with a set of predetermined realistic options to make sure it more or less works "
                      + "(default False). Ignores all other options; output file names can still be provided, "
                      + "otherwise the output files will start with a test_ prefix.")

    parser.add_option('-n','--number_of_samples', type='int', metavar='N', help="Number of samples to pool (required).")
    parser.add_option('-N','--number_of_pools', type='int', metavar='M', help="Number of resulting pools (required).")
    parser.add_option('-c','--binary_code_list_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a list of codewords).")
    parser.add_option('-C','--binary_code_generator_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a generator matrix).")

    parser.add_option('-m','--multiple_Biomek_files', action='store_true', default=True,
                      help="Generate multiple Biomek files, one per source plate (on by default).")
    parser.add_option('-o','--one_Biomek_file', action='store_false', dest='multiple_Biomek_files',
                      help="Generate a single Biomek file regardless of number of plates involved (off by default).")

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

    parser.add_option('-v','--volume_per_transfer', type='int', default=20, metavar='V', 
                      help="Liquid volume to use for each sample-to-pool transfer (default %default).")
    parser.add_option('-H','--Biomek_file_header', default="SourcePlt,SourceWell,DestPlt,DestWell,Volume", metavar='S', 
                      help="Header line for the Biomek transfer file (won't be printed if set to ''; default %default).")

    return parser


def get_binary_code(listfile=None,matrixfile=None):
    """ Given a listfile or matrixfile name (but not both!), return the generated binary code."""
    import binary_code_utilities
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

    # MAYBE-TODO could allow -p/-P to be automatically calculated from -n/-N and -s/-S?

    return options,outfile_basename


def generate_outfile_names(outfile_basename, if_multiple_files):
    """ Given the base outfile name, generate full outfile names (X -> X.txt and X_Biomek.csv or similar). """
    (outfile_data,outfile_Biomek) = [outfile_basename+suffix for suffix in ['.txt','_Biomek.csv']]
    # TODO implement the -m/-o options here - this should return a single outfile_Biomek with -o, or a list with -m.
    # TODO add this to the unit-test!
    return (outfile_data,outfile_Biomek)


def split_command_list_by_source(transfer_file_command_list):
    """ Given a list of "Source,_,_,_,_" strings, split it into multiple lists, one for each Source value."""
    # TODO implement
    # TODO add this to the unit-test!
    pass


def print_data_to_Biomek_files(transfer_file_command_list, outfile_Biomek, Biomek_header=""):
    ### Biomek file transfer file: just print the header and commands, nothing else
    save_line_list_as_file(transfer_file_command_list, outfile_Biomek, header=Biomek_header)
    # TODO need to come up with another mechanism with the -m option!


def print_data_to_outfile(sample_number_position_codeword_list, pool_number_position_list, outfile_data, options=None):
    """ Print data to main outfile: header information (command, path, date/time, options - all as #-start comments), 
    the sample_number_position_codeword_list, and the pool_number_position_list. """
    # print all the usual header information (command, path, date/time, options)
    OUTFILE = open(outfile_data,'w')
    write_header_data(OUTFILE,options)
    OUTFILE.write("# Corresponding Biomek file: %s\n"%outfile_Biomek)
    # write the sample_number_position_codeword_list and pool_number_position_list data
    OUTFILE.write("\nsample_number,\tplate_and_well_position,\tcodeword\n")
    for (n,p,c) in sample_number_position_codeword_list:
        OUTFILE.write("%s,\t%s,\t%s\n"%(n,p,c))
    OUTFILE.write("\npool_number,\tplate_and_well_position\n")
    for (n,p) in pool_number_position_list:
        OUTFILE.write("%s,\t%s\n"%(n,p))
    OUTFILE.close()


if __name__=='__main__':

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # Unit-testing - don't even look for more options/arguments, just run the test suite
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in simple test suite.")
        print "Defined plate sizes: %s"%plate_sizes
        test_functionality()
        # if not doing a test run, exit now; doing a normal run after unit-testing isn't allowed, but doing a test run is.
        if not options.test_run:
            sys.exit(0)

    # Test run: robotic_plate_transfer.py -n 63 -N 15 -P 3 -i Source1 -C error-correcting_codes/15-6-6_generator test
    # since we'll be redoing the parsing on an example string, the option value will be overwritten, so save it separately
    test_run = options.test_run     
    if test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test run.")
        # MAYBE-TODO could define multiple test runs later?
        test_input = "-n 63 -N 15 -P 3 -i Source1 -C error-correcting_codes/15-6-6_generator"
        # MAYBE-TODO the -C option value above will only work if we're in the directory where the script is - fix that?
        print "Test run arguments:", test_input
        # regenerate options with test argument string
        (options, _) = parser.parse_args(test_input.split())
        # if an outfile name was provided, accept it, otherwise set to 'test'
        if not args:    args = ['test']

    # run samples_and_code_to_Biomek_file
    options,outfile_basename = check_options_and_args(parser,options,args)
    all_outfilenames = generate_outfile_names(outfile_basename,options.multiple_Biomek_files)
    print "Output files:", all_outfilenames
    binary_code = get_binary_code(options.binary_code_list_file, options.binary_code_generator_file)
    args = (options.number_of_samples, options.number_of_pools, binary_code, options.volume_per_transfer, 
            options.number_of_sample_plates, options.size_of_sample_plates, options.sample_plate_IDs, 
            options.number_of_pool_plates, options.size_of_pool_plates, options.pool_plate_IDs)
    result_data_tuple = samples_and_code_to_Biomek_file(*args)
    (transfer_file_command_list, sample_number_position_codeword_list, pool_number_position_list) = result_data_tuple
    (outfile_data, outfile_Biomek) = all_outfilenames
    print_data_to_Biomek_files(transfer_file_command_list, outfile_Biomek, options.Biomek_file_header)
    print_data_to_outfile(sample_number_position_codeword_list, pool_number_position_list, outfile_data, options)

    if test_run:
        print("*** Test run finished. If you didn't get any errors, that's good. Check the output files to make sure.")

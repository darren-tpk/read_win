### Read WIN-format data (without ObsPy)
# Darren Tan, Gilles Seropian (March 2026)
# Note that ObsPy's read function does not handle WIN data appropriately. 
# Three main issues with ObsPy (for loading WIN data) have been reported by S. Nakagawa and A. Kato, namely (see https://www.eri.u-tokyo.ac.jp/people/nakagawa/win/):
# (1) bugs with data encoded 0.5- and 3-byte data
# (2) cannot deal with missing data or varying sampling rates
# (3) it is slow
# 
# These authors have provided a solution using the “libwinsystem.so”. Their solution seems to work well but requires the WIN system to be installed.
# Here we provide some code for people who cannot install the WIN system.

# References:
# https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/include/win/data_format.html

import time
import os
import glob
import warnings
import numpy as np
import pandas as pd
import struct
import obspy as obs
from obspy.core.trace import Stats

from pathlib import Path

class OneSecUnit(object):
    """
    Reads a 1 second unit into a Stream (one Trace per channel).
    
    Each 1 second unit is composed of:
        (1) A header containing the whole unit size and start time
        (2) A succession of channel blocks, each containing 1s of data for a given channel
    
    Attributes:
        cursor_position (int): Position of the reading cursor within the unit (in bytes)
        size (int): Size of the whole unit in bytes
        starttime (obs.UTCDateTime): Start time of the unit
        stream (obs.Stream): Stream containing data from unit
    """

    def __init__(self, unit_bytes):
        """
            Initializes class and automatically reads the 1 s unit.

        Args:
            unit_bytes (array of bytes): Bytes as read from WIN file. Assumes that 1s unit starts at position 0.
        """

        # Read header data
        self._read_header(unit_bytes)

        # Initial cursor reading position. Header always takes 10 bytes, then start of first channel block
        self.cursor_position = 10 

        # Initialize streams which will be populated with Traces as they are read
        self.stream = obs.Stream() 

        # Keep reading channel blocks until the end of the unit
        while self.cursor_position < self.size:
            self._read_channel_block(unit_bytes[self.cursor_position:self.size])


    def _read_header(self, unit_bytes):
        """Reading overall header of 1 s unit. Header contains two parts:
            (1) The first 4 bytes give the overall unit size in bytes (including these 4 bytes). Numbers >1byte are always big endian (most significant first).
            (2) The following 6 bytes give the start time of the unit, with each single byte giving year, month, day, hour, minute and second, respectively. For each single byte, the first 4 bits give the 10s digit andthe last 4 bits the 1s digit (i.e. 20 would be encoded as 0010-0000, that is 2-0, rather than 0010100, that is 20).
        
        Args:
            unit_bytes (array of bytes): Bytes as read directly from the winfile. 
        """
        self.size = struct.unpack(">I", unit_bytes[:4])[0] # First 4 bytes = total unit size in bytes
        self.starttime = _get_starttime(unit_bytes[4:10]) # Next 6 bytes are start time

    def _read_channel_block(self, channel_bytes):
        """Summary
        
        Args:
            channel_bytes (TYPE): Description
        """

        # First byte of channel block is channel ID
        ch_id = struct.unpack(">H", channel_bytes[:2])[0] 

        # Next two bytes are slightly tricky, first 0.5 byte is sample size in bytes then next 1.5 bytes is sampling_rate
        datasize_fs = struct.unpack(">H", channel_bytes[2:4])[0]

        datasize_code = datasize_fs >> 12 # Next 0.5 byte (i.e. 4 bits) is size of each data sample in bytes. If 0, then datasize is actually 0.5 byte
        sampling_rate = datasize_fs & 4095 # Next 1.5 bytes (i.e. 12 bits) is sampling rate in Hz. Note that 4095 = 0xFFF

        # If datasize is 0, then actual datasize is 0.5 byte
        if datasize_code == 0:
            datasize = 0.5
        else:
            datasize = datasize_code

        # Determine SEED channel code based on sampling rate
        if 10 <= sampling_rate < 80:
            channel_code = 'BDF'
        elif 80 <= sampling_rate < 250:
            channel_code = 'HDF'
        elif 250 <= sampling_rate < 1000:
            channel_code = 'CDF'
        else:
            warnings.warn('Sampling rate is < 10 or >= 1000 Hz!')

        # Prepare metadata for trace (note that more metadata, e.g. network, will be added later)
        stats = Stats({
                        "sampling_rate": sampling_rate,
                        "npts": sampling_rate,
                        "location": f'{ch_id:04d}',
                        "channel": channel_code, 
                        "starttime": self.starttime,
                    })

        # First data sample, always 4 bytes, in absolute units 
        y_0 = uint2int( struct.unpack(">I", channel_bytes[4:8])[0], 4)

        # Format strings controls how bytes are unpacked. Correct letter will be selected based on datasize
        format_string = 'BBHBI'

        # Read all remaining bytes at once, but for 0.5 and 3 byte sample size, things are slightly trickier.

        if datasize_code == 0: # for 0.5 byte size, data samples are 4-bit but we read 8-bit, so each byte read actually contains two data samples
            N_samples_to_read = (sampling_rate - 1) // 2 
            N_bytes_to_read = (sampling_rate - 1) // 2   

            # If sampling rate is even, e.g. 100 Hz, there are an odd numbers of data samples to read (e.g. 99). But the data size is 0.5 bytes, so the overall size will not be an integer (e.g. 49.5 bytes)
            # In this case, there will be an additional empty meaningless "padding" 0.5 bytes (i.e. 4 bits) at the end of the 
            if  (sampling_rate - 1) % 2:
                N_samples_to_read += 1
                N_bytes_to_read += 1

            diff_udata_8bit = np.array( struct.unpack('>' + format_string[datasize_code] * N_samples_to_read, channel_bytes[8:8+N_bytes_to_read]) ) 
            
            first_4bit = diff_udata_8bit >> 4
            second_4bit = diff_udata_8bit & 0xf
            diff_udata = np.stack( (first_4bit, second_4bit), axis=1).flatten()

            # For even sampling rates, the last 4 bits are zeros for padding only and not actual data
            if  (sampling_rate - 1) % 2:
                diff_udata = diff_udata[:-1]

        elif datasize_code == 3: # for 3 byte size, data samples are 24-bit which cannot be read. Instead, we read 8-bit, and combine 3 bytes to get a data sample.
            N_samples_to_read = (sampling_rate - 1) * 3 
            N_bytes_to_read = N_samples_to_read 
            diff_udata_8bit = np.array( struct.unpack('>' + format_string[datasize_code] * N_samples_to_read, channel_bytes[8:8+N_bytes_to_read]) ) 

            # Recombine each set of three 8-bit value into single 24-bit values
            first_8bit, second_8bit, third_8bit = np.reshape(diff_udata_8bit, (diff_udata_8bit.shape[0] // 3, 3)).T
            diff_udata = (first_8bit << 16) + (second_8bit << 8) + third_8bit

        else:
            N_samples_to_read = sampling_rate - 1 # number of data samples remaining to be read
            N_bytes_to_read = N_samples_to_read * datasize # number of actual bytes to be read
            diff_udata = np.array( struct.unpack('>' + format_string[datasize_code] * N_samples_to_read, channel_bytes[8:8+N_bytes_to_read]) ) # reading all remaining data at once. Note that these are differential, unsigned integers

        # Convert from unsigned to signed integers, still differential (i.e. difference to previous value rather than absolute value)
        diff_data = to_signed_array(diff_udata, datasize)

        # Convert to absolute values and add to unit Stream 
        # Note that the data is still in integers, and not in physical units yet (will be converted all at once, after all 1s units are loaded)
        self.stream += obs.Trace(data=np.cumsum( np.insert(diff_data, 0, y_0) ).astype(float), header=stats)

        # Update cursor position to end of this channel block
        self.cursor_position += int( 8 + (sampling_rate - 1) * datasize )

        # If (sampling_rate - 1) is odd, the previous line of code will round down and an extra byte needs to be added for padding
        if ( datasize_code == 0 ) & ( (sampling_rate - 1) % 2 ):
            self.cursor_position += 1


def read_header(all_bytes):
    """
    Reading overall header of 1 s unit. Header contains two parts:
        (1) The first 4 bytes give the overall unit size in bytes (including these 4 bytes). Numbers >1byte are always big endian (most significant first).
        (2) The following 6 bytes give the start time of the unit, with each single byte giving year, month, day, hour, minute and second, respectively. For each single byte, the first 4 bits give the 10s digit andthe last 4 bits the 1s digit (i.e. 20 would be encoded as 0010-0000, that is 2-0, rather than 0010100, that is 20).
    
    Args:
        all_bytes (array of bytes): Bytes as read directly from the winfile. 
    """
    unit_size = struct.unpack(">I", all_bytes[:4])[0] # First 4 bytes = total unit size in bytes
    starttime = _get_starttime(all_bytes[4:10]) # Next 6 bytes are start time

    return unit_size, starttime

def _get_starttime(date_bytes):
    """Extracts the start time from a 1 s unit.
    
    Args:
        date_bytes (array): The 6 bytes giving the date in the header.
    
    Returns:
        starttime (obs.UTCDateTime): Start time of 1 s unit as obs.UTCDateTime
    """
    starttime = [_dateint2num(i) for i in date_bytes] # Careful, when using list comprehension, each item b from all_bytes will actually be an integer rather than a byte

    if starttime[0] <= 80:
        starttime[0] += 2000 # add 2000 to years between 2000 and 2080
    else:
        starttime[0] += 1900 # add 1900 to years between 1981 and 1999 

    return obs.UTCDateTime(*starttime)

def _dateint2num(date_int):
    # Takes an integer and decodes it into its date value
    # First 4 bits are for the decade digit and last 4 bits are for the unit digit
    return ( date_int >>  4 ) * 10 + ( date_int &  15 ) # Note that 15 = 0xF

def _datebyte2num(byte):
    # Takes a byte and decodes it into its date value
    # First 4 bits are for the decade digit and last 4 bits are for the unit digit
    return (struct.unpack("B", byte)[0] >>  4) * 10 + (struct.unpack("B", byte)[0] &  15)

def read_channel_table(channel_table_path):
    # Read channel table into a pandas DataFrame
    if not Path(channel_table_path).exists():
        raise FileNotFoundError("Channel table path does not exist.")

    print(f'Loading channel table from {channel_table_path}...', end='')

    COLUMNS = columns = ["location", "flag", "delay_time", "station", "channel",
                         "amplitude_reduction", "digitizing_bits",  "sensitivity",
                         "unit", "natural_period", "damping_constant", "amplification_factor",
                         "digitization_magnitude", "latitude", "longitude", "altitude",
                         "correction_P", "correction_S"]
    channel_table = pd.read_csv(channel_table_path, sep=r"\s+", header=None, dtype={0: str})
    channel_table.columns = columns[:channel_table.shape[1]]

    # Calculate amplitude correction factor
    channel_table["amplitude_correction"] = (channel_table["digitization_magnitude"] *
                                            (10 ** (channel_table["amplification_factor"] / 20)) /
                                            channel_table["sensitivity"])

    print('   Done!')

    return channel_table

def to_signed_array(uarray, byte_size):
    """Converts an array of unsigned integers to signed integers.

    For instance, if the byte size is 2, the number of bits is 16. Unsigned integers will range from 0 to 2^16 - 1 = 65535. The corresponding signed integers will range from -32767 to 32768.
    
    Args:
        uarray (array): Array of unsigned integer values
        byte_size (int): Size of integer in bytes. E.g. byte_size = 2 means 16-bit integer.
    
    Returns:
        sarray (array): Array of corresponding signed integer values
    """
    # Number of bits
    n_bits = int( byte_size * 8 )

    # Check unsigned is in the correct range
    if np.any( (uarray < 0) | ( uarray >= 2 ** (n_bits)) ):
        raise ValueError(f'Unsigned array contains values outside of range for {n_bits}-bit integers. Values should be between 0 and {2 ** (n_bits) - 1}.')

    # Convert whole array to signed integers
    sarray = uarray.copy()
    sarray[ sarray > 2 ** ( n_bits - 1 )] -= 2 ** n_bits

    return sarray

def uint2int(unsigned, byte_size):
    """Given an unsigned integer value (from 0 to 2^n_bits - 1), return the corresponding signed integer (from -2^(n_bits-1)+1 to 2^(n_bits-1).

    For instance, if the byte size is 2, the number of bits is 16. Unsigned integers will range from 0 to 2^16 - 1 = 65535. The corresponding signed integers will range from -32767 to 32768.
    
    Args:
        unsigned (int): Unsigned integer value
        byte_size (int): Size of integer in bytes. E.g. byte_size = 2 means 16-bit integer.
    
    Returns:
        signed (int): Unsigned integer value
    """
    # Number of bits
    n_bits = byte_size * 8

    # Check unsigned is in the correct range
    if ( unsigned < 0 ) | ( unsigned >= 2 ** (n_bits)):
        raise ValueError(f'Unsigned integer {unsigned} outside of range for {n_bits}-bit integers. Values should be between 0 and {2 ** (n_bits) - 1}.')

    elif unsigned <= 2 ** ( n_bits - 1 ):
        signed = unsigned

    else:
        signed = unsigned - 2 ** n_bits 

    return signed

def get_opposite(number, byte_size):
    """Return the 2-complementary number, for a given byte_size
    
    Args:
        number (int): Integer to be "inverted"
        byte_size (int): Byte size
    
    Returns:
        complementary (int): 2-complementary of number for given byte_size. E.g. for 16-bit integers, if number = 1234 (0100 1101 0010), complementary = 64302 (1111 1011 0010 1110).
    """

    # Number of bits
    n_bits = byte_size * 8
    complementary = 2 ** n_bits - number

    return complementary

def _read_single_win_file(file_path):
    """Reads a single WIN file into an Obspy.Stream object.
    
    Args:
        file_path (str): Path to WIN file
    
    Returns:
        stream (obs.Stream): Obspy Stream object containing data from WIN file (not merged or trimmed)
    """

    # Print statement to let user know something is happening!
    print(f'.......................Reading WIN file {file_path}...', end='')

    # Read all bytes from winfile into a variable
    with open(file_path, "rb") as f:
        all_bytes = f.read()

    # Initialize Stream to load data into
    stream = obs.Stream()

    # Position of cursor as reading through all 1 s units
    read_position = 0

    # Keep reading 1s units until the end of the file
    while read_position < len(all_bytes):

        unit = OneSecUnit( all_bytes[read_position:] ) # read one second unit

        stream += unit.stream # add data to stream

        read_position += unit.size # advance cursor position to start of next unit

    print(' Done !')

    return stream


def read_win(list_file_paths, channel_table_path, network='SA', station='KUR', fill_value=None):
    """Reads WIN file(s) into an Obspy Stream.
    
    Args:
        list_file_paths (str or list): Path or list of paths to WIN file
        channel_table_path (str): Path to channel table
        network (str, optional): SEED network code (2 characters, A-Z)
        station (str, optional): SEED station code (3-4 characters, A-Z & 0-9)
        fill_value (optional): Fill value for merging Traces with gaps
    
    Returns:
        stream (obs.Stream): Stream containing infrasound data
    """

    run_clock = time.time()

    if len(list_file_paths) == 0:
        raise FileNotFoundError('Input list of files is empty.')
    
    elif type(list_file_paths) is str: # if single path given as str, turn it into a list
            list_file_paths = [list_file_paths]

    print(f'Beginning conversion of {len(list_file_paths)} file(s) to Stream...')

    # Initialize Stream to load data into
    stream = obs.Stream()

    # Load each file in the list one by one
    for win_file in list_file_paths:
        stream += _read_single_win_file(win_file)
        # stream.merge(fill_value=fill_value) # could merge after each iteration or once at the end... Not sure which one is faster or if it makes a difference.

    # print('All files successfully read and streams merged. Applying amplitude correction from channel table.')

    # Could also merge only once here... Seems slightly faster
    print('-----All files successfully read. Merging all streams...')
    stream.merge(fill_value=fill_value)
    print('-----Streams merged. Applying amplitude correction from channel table.')


    # Load channel table
    channel_table = read_channel_table(channel_table_path)

    # Correct amplitude with values from channel table and add metadata
    for tr in stream:
        tr.stats.network = network
        tr.stats.station = station
        tr.data *= channel_table.loc[ channel_table['location'] == str(tr.stats.location), 'amplitude_correction' ].values[0]

    print('...finished conversion process to Stream object! Time elapsed: %s seconds.\n' % (time.time() - run_clock))

    return stream


if __name__ == '__main__':

    WIN_FILE_PATH = './kurokami123/05/2603051*' # either a single file path or can glob to a pattern
    CHANNEL_TABLE_PATH = './channels.tbl'

    list_file_paths = glob.glob(WIN_FILE_PATH)
    list_file_paths.sort() # sorting files not actually needed

    st = read_win(list_file_paths, CHANNEL_TABLE_PATH, network='SA', station='KUR', fill_value=None)
import os
import warnings
import numpy as np
import pandas as pd
import struct
import warnings
from glob import glob
from obspy import UTCDateTime, Stream, Trace
from obspy.core.trace import Stats
from pathlib import Path

class OneSecUnit(object):
    """
    Read a 1-second WIN data unit into an ObsPy Stream (one Trace per channel).

    Each unit consists of:
    (1) A header containing the unit size and start time
    (2) A sequence of channel blocks, each containing 1 second of data

    :param unit_bytes (bytes or array): Raw bytes for a single 1-second unit

    :ivar cursor_position (int): Current reading position within the unit (in bytes)
    :ivar size (int): Total size of the unit in bytes
    :ivar starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time of the unit
    :ivar stream (:class:`~obspy.core.stream.Stream`): Stream containing data from the unit
    """

    def __init__(self, unit_bytes):
        """
        Initialize and parse a 1-second WIN data unit.

        :param unit_bytes (bytes or array): Raw bytes for a single 1-second unit
        :return: None
        """

        # Read header data
        self._read_header(unit_bytes)

        # Initial cursor reading position. Header always takes 10 bytes, then start of first channel block
        self.cursor_position = 10

        # Initialize streams which will be populated with Traces as they are read
        self.stream = Stream()

        # Keep reading channel blocks until the end of the unit
        while self.cursor_position < self.size:
            self._read_channel_block(unit_bytes[self.cursor_position:self.size])

    def _read_header(self, unit_bytes):
        """
        Read the header of a 1-second WIN data unit.
        The header consists of:
        (1) 4-byte big-endian integer giving total unit size (including header)
        (2) 6 bytes encoding start time (YY MM DD HH MM SS), where each byte
            stores two decimal digits (upper 4 bits = tens, lower 4 bits = units)

        :param unit_bytes (bytes or array): Raw bytes from the WIN file
        :return: None
        """
        self.size = struct.unpack(">I", unit_bytes[:4])[0]  # First 4 bytes = total unit size in bytes
        self.starttime = _get_starttime(unit_bytes[4:10])  # Next 6 bytes are start time

    def _read_channel_block(self, channel_bytes):
        """
        Read a single channel block and append it as a Trace to the Stream.
        Handles special cases for 0.5-byte (4-bit) and 3-byte (24-bit) data encoding.
        Each channel block contains:
        - Channel ID
        - Sample size and sampling rate (bit-packed)
        - One absolute sample followed by differential samples

        :param channel_bytes (bytes or array): Bytes corresponding to a single channel block
        :return: None
        """

        # First byte of channel block is channel ID
        ch_id = struct.unpack(">H", channel_bytes[:2])[0]

        # Next two bytes are slightly tricky, first 0.5 byte is sample size in bytes then next 1.5 bytes is sampling_rate
        datasize_fs = struct.unpack(">H", channel_bytes[2:4])[0]

        datasize_code = datasize_fs >> 12  # Next 0.5 byte (i.e. 4 bits) is size of each data sample in bytes. If 0, then datasize is actually 0.5 byte
        sampling_rate = datasize_fs & 4095  # Next 1.5 bytes (i.e. 12 bits) is sampling rate in Hz. Note that 4095 = 0xFFF

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
        y_0 = uint2int(struct.unpack(">I", channel_bytes[4:8])[0], 4)

        # Format strings controls how bytes are unpacked. Correct letter will be selected based on datasize
        format_string = 'BBHBI'

        # Read all remaining bytes at once, but for 0.5 and 3 byte sample size, things are slightly trickier.

        if datasize_code == 0:  # for 0.5 byte size, data samples are 4-bit but we read 8-bit, so each byte read actually contains two data samples
            N_samples_to_read = (sampling_rate - 1) // 2
            N_bytes_to_read = (sampling_rate - 1) // 2

            # If sampling rate is even, e.g. 100 Hz, there are an odd numbers of data samples to read (e.g. 99). But the data size is 0.5 bytes, so the overall size will not be an integer (e.g. 49.5 bytes)
            # In this case, there will be an additional empty meaningless "padding" 0.5 bytes (i.e. 4 bits) at the end of the
            if (sampling_rate - 1) % 2:
                N_samples_to_read += 1
                N_bytes_to_read += 1

            diff_udata_8bit = np.array(struct.unpack('>' + format_string[datasize_code] * N_samples_to_read,
                                                     channel_bytes[8:8 + N_bytes_to_read]))

            first_4bit = diff_udata_8bit >> 4
            second_4bit = diff_udata_8bit & 0xf
            diff_udata = np.stack((first_4bit, second_4bit), axis=1).flatten()

            # For even sampling rates, the last 4 bits are zeros for padding only and not actual data
            if (sampling_rate - 1) % 2:
                diff_udata = diff_udata[:-1]

        elif datasize_code == 3:  # for 3 byte size, data samples are 24-bit which cannot be read. Instead, we read 8-bit, and combine 3 bytes to get a data sample.
            N_samples_to_read = (sampling_rate - 1) * 3
            N_bytes_to_read = N_samples_to_read
            diff_udata_8bit = np.array(struct.unpack('>' + format_string[datasize_code] * N_samples_to_read,
                                                     channel_bytes[8:8 + N_bytes_to_read]))

            # Recombine each set of three 8-bit value into single 24-bit values
            first_8bit, second_8bit, third_8bit = np.reshape(diff_udata_8bit, (diff_udata_8bit.shape[0] // 3, 3)).T
            diff_udata = (first_8bit << 16) + (second_8bit << 8) + third_8bit

        else:
            N_samples_to_read = sampling_rate - 1  # number of data samples remaining to be read
            N_bytes_to_read = N_samples_to_read * datasize  # number of actual bytes to be read
            diff_udata = np.array(struct.unpack('>' + format_string[datasize_code] * N_samples_to_read, channel_bytes[
                                                                                                        8:8 + N_bytes_to_read]))  # reading all remaining data at once. Note that these are differential, unsigned integers

        # Convert from unsigned to signed integers, still differential (i.e. difference to previous value rather than absolute value)
        diff_data = to_signed_array(diff_udata, datasize)

        # Convert to absolute values and add to unit Stream
        # Note that the data is still in integers, and not in physical units yet (will be converted all at once, after all 1s units are loaded)
        self.stream += Trace(data=np.cumsum(np.insert(diff_data, 0, y_0)).astype(float), header=stats)

        # Update cursor position to end of this channel block
        self.cursor_position += int(8 + (sampling_rate - 1) * datasize)

        # If (sampling_rate - 1) is odd, the previous line of code will round down and an extra byte needs to be added for padding
        if (datasize_code == 0) & ((sampling_rate - 1) % 2):
            self.cursor_position += 1


def _get_starttime(date_bytes):
    """
    Extract the start time from a 1-second unit header.

    :param date_bytes (array): The 6 bytes encoding the date in the header
    :return: :class:`~obspy.core.utcdatetime.UTCDateTime`
    """
    starttime = [_dateint2num(i) for i in
                 date_bytes]  # Careful, when using list comprehension, each item b from all_bytes will actually be an integer rather than a byte

    if starttime[0] <= 80:
        starttime[0] += 2000  # add 2000 to years between 2000 and 2080
    else:
        starttime[0] += 1900  # add 1900 to years between 1981 and 1999

    return UTCDateTime(*starttime)


def _dateint2num(date_int):
    """
    Decode an integer into its corresponding date value.
    The input integer is interpreted as two packed decimal digits:
    the upper 4 bits represent the decade digit, and the lower
    4 bits represent the unit digit (where 0xF = 15).

    :param date_int (int): Integer to be decoded
    :return: int
    """
    return (date_int >> 4) * 10 + (date_int & 15)


def _datebyte2num(byte):
    """
    Decode a byte into its corresponding date value.
    The input byte is interpreted as two packed decimal digits:
    the upper 4 bits represent the decade digit, and the lower
    4 bits represent the unit digit.

    :param byte (bytes): Single byte to be decoded
    :return: int
    """
    return (struct.unpack("B", byte)[0] >> 4) * 10 + (struct.unpack("B", byte)[0] & 15)


def read_channel_table(channel_table_path):
    """
    Read a channel table file into a pandas DataFrame.

    :param channel_table_path (str): Path to channel table file
    :return: :class:`pandas.DataFrame`
    """

    # Read channel table into a pandas DataFrame
    if not Path(channel_table_path).exists():
        raise FileNotFoundError("Channel table path does not exist.")

    print(f'Loading channel table from {channel_table_path}...', end='')

    COLUMNS = columns = ["location", "flag", "delay_time", "station", "channel",
                         "amplitude_reduction", "digitizing_bits", "sensitivity",
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
    """
    Convert an array of unsigned integers to signed integers for a given byte size.
    For example, for byte_size = 2 (16-bit), unsigned integers range from 0 to 65535,
    and the corresponding signed integers range from -32767 to 32768.

    :param uarray (array): Array of unsigned integer values
    :param byte_size (int): Size of integer in bytes (e.g., 2 for 16-bit)
    :return: array
    """
    # Number of bits
    n_bits = int(byte_size * 8)

    # Check unsigned is in the correct range
    if np.any((uarray < 0) | (uarray >= 2 ** (n_bits))):
        raise ValueError(
            f'Unsigned array contains values outside of range for {n_bits}-bit integers. Values should be between 0 and {2 ** (n_bits) - 1}.')

    # Convert whole array to signed integers
    sarray = uarray.copy()
    sarray[sarray > 2 ** (n_bits - 1)] -= 2 ** n_bits

    return sarray


def uint2int(unsigned, byte_size):
    """
    Convert an unsigned integer to its signed representation for a given byte size.

    Given an unsigned integer in the range [0, 2^n_bits - 1], return the corresponding
    signed integer. For example, for byte_size = 2 (16-bit), unsigned integers range
    from 0 to 65535, and the corresponding signed integers range from -32767 to 32768.

    :param unsigned (int): Unsigned integer value
    :param byte_size (int): Size of integer in bytes (e.g., 2 for 16-bit)
    :return: int
    """
    # Number of bits
    n_bits = byte_size * 8

    # Check unsigned is in the correct range
    if (unsigned < 0) | (unsigned >= 2 ** (n_bits)):
        raise ValueError(
            f'Unsigned integer {unsigned} outside of range for {n_bits}-bit integers. Values should be between 0 and {2 ** (n_bits) - 1}.')

    elif unsigned <= 2 ** (n_bits - 1):
        signed = unsigned

    else:
        signed = unsigned - 2 ** n_bits

    return signed


def get_opposite(number, byte_size):
    """
    Return the two's complement of an integer for a given byte size.
    For example, for 16-bit integers, if number = 1234 (0100 1101 0010),
    the complementary value is 64302 (1111 1011 0010 1110).

    :param number (int): Integer to be converted
    :param byte_size (int): Byte size of the integer representation
    :return: int
    """

    # Number of bits
    n_bits = byte_size * 8
    complementary = 2 ** n_bits - number

    return complementary


def _read_single_win_file(file_path):
    """
    Read a single WIN file into an ObsPy Stream object (not merged or trimmed)

    :param file_path (str): Path to WIN file
    :return: :class:`~obspy.core.stream.Stream`
    """

    # Print statement to let user know something is happening!
    print(f'.......................Reading WIN file {file_path}...', end='')

    # Read all bytes from winfile into a variable
    with open(file_path, "rb") as f:
        all_bytes = f.read()

    # Initialize Stream to load data into
    stream = Stream()

    # Position of cursor as reading through all 1 s units
    read_position = 0

    # Keep reading 1s units until the end of the file
    while read_position < len(all_bytes):
        unit = OneSecUnit(all_bytes[read_position:])  # read one second unit

        stream += unit.stream  # add data to stream

        read_position += unit.size  # advance cursor position to start of next unit

    print(' Done !')

    return stream

def _safe_merge(stream, fill_value):
    """
    Merge Traces with the same ID, modifying data types and rounding non-integer sampling rates if necessary.
    Modified from code by Aaron Wech; ported over from uafgeotool's waveform_collection repository.

    :param stream (:class:`~obspy.core.stream.Stream`): Input Stream (modified in-place)
    :param fill_value (int, float, str, or None): Passed to :meth:`obspy.core.stream.Stream.merge`
    :return: None
    """

    try:
        stream.merge(fill_value=fill_value)
    except Exception:  # ObsPy raises an Exception if data types are not all identical
        for trace in stream:
            if trace.data.dtype != np.dtype(np.int32):
                trace.data = trace.data.astype(np.int32, copy=False)
    try:
        stream.merge(fill_value=fill_value)
    except Exception:  # ObsPy also raises an Exception if traces with the same ids have different sampling rates
        for trace in stream:
            if trace.stats.sampling_rate != np.round(trace.stats.sampling_rate):
                warnings.warn('Rounding off %s sampling rate from %f Hz to %.1f Hz for merge compatibility.' % (
                              trace.id, trace.stats.sampling_rate, np.round(trace.stats.sampling_rate)))
                trace.stats.sampling_rate = np.round(trace.stats.sampling_rate)
        stream.merge(fill_value=fill_value)  # Try merging with rounded sampling rates

def read_win(file_paths, channel_table_path, fill_value=None):
    """
    Read WIN file(s) into an ObsPy Stream.

    :param file_paths (str or list): Path or list of paths to WIN file(s)
    :param channel_table_path (str): Path to channel table file
    :param fill_value (float, optional): Fill value used in stream.merge() to handle gaps
    :return: :class:`~obspy.core.stream.Stream`
    """

    # Check input file paths
    if len(file_paths) == 0:
        raise FileNotFoundError('Input list of files is empty.')
    elif type(file_paths) is str:  # if single path given as str, turn it into a list
        file_paths = [file_paths]
    print(f'Beginning conversion of {len(file_paths)} file(s) to Stream...')

    # Initialize Stream to load data into
    stream = Stream()

    # Load each file in the list one by one
    for win_file in file_paths:
        stream += _read_single_win_file(win_file)

    # print('All files successfully read and streams merged. Applying amplitude correction from channel table.')

    # Could also merge only once here... Seems slightly faster
    print('-----All files successfully read. Merging all streams...')
    _safe_merge(stream, fill_value)
    print('-----Streams merged. Applying amplitude correction from channel table.')

    # Load channel table
    channel_table = read_channel_table(channel_table_path)

    # Correct amplitude with values from channel table and add metadata
    for trace in stream:
        trace.stats.network = network
        trace.stats.station = station
        trace.data *= \
        channel_table.loc[channel_table['location'] == str(trace.stats.location), 'amplitude_correction'].values[0]

    print('...finished conversion process to Stream object!\n')

    return stream


def read_win_dir(starttime,
                 endtime,
                 file_directory,
                 file_pattern,
                 file_interval,
                 fill_value,
                 channel_table_path):
    """
    Read WIN format data from a directory or archive.

    :param starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time of desired data
    :param endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): End time of desired data
    :param file_directory (str): Path to directory containing WIN data
    :param file_pattern (str): WIN data file pattern (must be a valid strftime input for UTCDateTime)
    :param file_interval (str): WIN data file interval ("minute", "hour", or "day")
    :param fill_value (float): Value passed to stream.merge() when combining data from separate WIN files
    :param channel_table_path (str): Path to channel table file
    :return: :class:`~obspy.core.stream.Stream`
    """

    # Check if file interval is a valid input:
    if file_interval not in ["minute", "hour", "day"]:
        raise ValueError("Invalid value for file_interval -- please use 'minute' or 'hour' or 'day'")

    # Get comprehensive list of WIN file times to be read
    if file_interval == "minute":
        file_interval_time = 60
        first_file_time = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour, starttime.minute)
        end_file_time = UTCDateTime(endtime.year, endtime.month, endtime.day, endtime.hour, endtime.minute)
    elif file_interval == "hour":
        file_interval_time = 3600
        first_file_time = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour)
        end_file_time = UTCDateTime(endtime.year, endtime.month, endtime.day, endtime.hour)
    else:
        file_interval_time = 86400
        first_file_time = UTCDateTime(starttime.year, starttime.month, starttime.day)
        end_file_time = UTCDateTime(endtime.year, endtime.month, endtime.day)
    if end_file_time == endtime:
        file_times = np.arange(first_file_time, end_file_time + file_interval_time, file_interval_time)
    else:
        file_times = np.arange(first_file_time, end_file_time + file_interval_time, file_interval_time)

    # Get list of WIN file paths
    file_paths = []
    for file_time in file_times:
        pattern = os.path.join(file_directory, file_time.strftime(file_pattern))
        matching_files = glob(pattern)
        file_paths.extend(matching_files)

    # Read WIN and trim
    stream = read_win(file_paths, channel_table_path, fill_value=fill_value)
    stream = stream.trim(starttime, endtime)

    return stream


if __name__ == '__main__':

    starttime = UTCDateTime(2026, 3, 5, 4, 30)  # UTC
    endtime = UTCDateTime(2026, 3, 5, 14, 15)  # UTC
    file_directory = "./kurokami*/*/"  # for glob
    file_pattern = "%y%m%d%H"  # UTCDateTime strftime input
    file_interval = "hour"  # "minute", "hour", or "day"
    fill_value = None  # for stream.merge
    channel_table_path = "./channels.tbl"

    stream = read_win_dir(starttime,
                          endtime,
                          file_directory,
                          file_pattern,
                          file_interval,
                          fill_value,
                          channel_table_path)
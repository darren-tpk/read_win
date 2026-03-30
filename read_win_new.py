import os
import time
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
    :ivar size_bytes (int): Total size of the unit in bytes
    :ivar starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time of the unit
    :ivar data (:class:`numpy.ndarray`): Array containing decoded waveform samples for the unit
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

        # Initialize placeholder for data
        self.data = {}

        # Keep reading channel blocks until the end of the unit
        while self.cursor_position < self.size_bytes:
            self._read_channel_block(unit_bytes[self.cursor_position:self.size_bytes])

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
        self.size_bytes = struct.unpack(">I", unit_bytes[:4])[0]  # First 4 bytes = total unit size in bytes
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

        # First byte of channel block is sensor ID
        sensor_id = f"{struct.unpack('>H', channel_bytes[:2])[0]:04x}"

        # Next two bytes are slightly tricky, first 0.5 byte is sample size in bytes then next 1.5 bytes is sampling_rate
        data_properties_bytes = struct.unpack(">H", channel_bytes[2:4])[0]
        datasize_code = data_properties_bytes >> 12  # Next 0.5 byte (i.e. 4 bits) is size of each data sample in bytes. If 0, then datasize is actually 0.5 byte
        datasize = 0.5 if datasize_code == 0 else datasize_code # If datasize is 0, then actual datasize is 0.5 byte
        sampling_rate = data_properties_bytes & 4095  # Next 1.5 bytes (i.e. 12 bits) is sampling rate in Hz. Note that 4095 = 0xFFF

        # Initialize data entry if not present
        if sensor_id not in self.data:
            self.data[sensor_id] = {"blocks": []}

        # First data sample, always 4 bytes, in absolute units
        first_sample_value = uint2int(struct.unpack(">I", channel_bytes[4:8])[0], 4)

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
            diff_udata = np.array(struct.unpack('>' + format_string[datasize_code] * N_samples_to_read, channel_bytes[8:8 + N_bytes_to_read]))  # reading all remaining data at once. Note that these are differential, unsigned integers

        # Convert from unsigned to signed integers, still differential (i.e. difference to previous value rather than absolute value)
        diff_data = to_signed_array(diff_udata, datasize)

        # Convert to absolute units
        data_segment = np.cumsum(np.insert(diff_data, 0, first_sample_value)).astype(float)

        # Append data segment
        # Note that the data is still in integers, and not in physical units yet (will be converted all at once, after all 1s units are loaded)
        self.data[sensor_id]["blocks"].append({"starttime": self.starttime,"sampling_rate": sampling_rate, "data": data_segment})

        # Update cursor position to end of this channel block
        self.cursor_position += int(8 + (sampling_rate - 1) * datasize)

        # If (sampling_rate - 1) is odd, the previous line of code will round down and an extra byte needs to be added for padding
        if (datasize_code == 0) & ((sampling_rate - 1) % 2):
            self.cursor_position += 1


def _get_starttime(date_bytes):
    """
    Extract and decode the start time from a 1-second unit header.

    Each input integer represents a time component encoded as two packed
    decimal digits, where the upper 4 bits correspond to the tens digit
    and the lower 4 bits correspond to the units digit (where 0xF = 15).

    :param date_bytes (array): The 6 bytes encoding the date in the header
    :return: :class:`~obspy.core.utcdatetime.UTCDateTime`
    """
    starttime = [(b >> 4) * 10 + (b & 0x0F) for b in date_bytes]  # Careful, when using list comprehension, each item b from all_bytes will actually be an integer rather than a byte
    starttime[0] += 2000 if starttime[0] <= 80 else 1900

    return UTCDateTime(*starttime)


def read_channel_table(channel_table_path):
    """
    Read a channel table file into a pandas DataFrame.

    :param channel_table_path (str): Path to channel table file
    :return: :class:`pandas.DataFrame`
    """

    # Read channel table into a pandas DataFrame
    if not Path(channel_table_path).exists():
        raise FileNotFoundError("Channel table path does not exist.")

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


def _read_single_win_file(file_path, verbose=True):
    """
    Read a single WIN file into an ObsPy Stream object (not merged or trimmed)

    :param file_path (str): Path to WIN file
    :param verbose (bool): If `False`, all print statements will be blocked. Default is `True`.
    :return: :class:`~obspy.core.stream.Stream`
    """

    # Use print for logging if verbose; otherwise use a no-op function
    log = print if verbose else lambda *args, **kwargs: None

    # Print statement to let user know something is happening!
    log(f'.......................Reading WIN file {file_path}...', end='')

    # Read all bytes from winfile into a variable
    with open(file_path, "rb") as f:
        all_bytes = f.read()

    # Initialize dictionary to accumulate data across all 1 s units
    all_data = {}

    # Position of cursor as reading through all 1 s units
    read_position = 0

    # Keep reading 1 s units until the end of the file
    while read_position < len(all_bytes):
        unit = OneSecUnit(all_bytes[read_position:])

        for sensor_id, sensor_info in unit.data.items():
            if sensor_id not in all_data:
                all_data[sensor_id] = {"blocks": []}

            all_data[sensor_id]["blocks"].extend(sensor_info["blocks"])

        read_position += unit.size_bytes

    # Initialize assembled data for WIN file
    assembled_data = {}

    for sensor_id, sensor_dict in all_data.items():
        blocks = sensor_dict["blocks"]

        # Sort blocks by starttime
        blocks = sorted(blocks, key=lambda b: b["starttime"])

        # Determine expected sampling rate (mode)
        sampling_rates = [b["sampling_rate"] for b in blocks]
        expected_sampling_rate = max(set(sampling_rates), key=sampling_rates.count)

        # Determine master start and end times
        master_starttime = blocks[0]["starttime"]
        master_endtime = max(b["starttime"] + len(b["data"]) / expected_sampling_rate for b in blocks)

        # Allocate master array
        total_samples = int(round((master_endtime - master_starttime) * expected_sampling_rate))
        master_data = np.full(total_samples, np.nan, dtype=float)

        # Fill data
        for b in blocks:

            # Check sampling rate
            if b["sampling_rate"] != expected_sampling_rate:
                warnings.warn(f"Skipping block for sensor {sensor_id} at {b['starttime']} "f"due to sampling rate mismatch ({b['sampling_rate']} != {expected_sampling_rate})")
                continue

            # Check number of samples
            if len(b["data"]) != expected_sampling_rate:
                warnings.warn(f"Skipping block for sensor {sensor_id} at {b['starttime']} "f"due to unexpected sample count ({len(b['data'])} != {expected_sampling_rate})")
                continue

            # Determine indices and check for overlap
            start_index = int(round((b["starttime"] - master_starttime) * expected_sampling_rate))
            end_index = start_index + expected_sampling_rate
            existing = master_data[start_index:end_index]
            if np.any(~np.isnan(existing)):
                warnings.warn(f"Overlap detected for sensor {sensor_id} at {b['starttime']} — existing data will be overwritten")

            # Write data
            master_data[start_index:end_index] = b["data"]

        assembled_data[sensor_id] = {"starttime": master_starttime, "sampling_rate": expected_sampling_rate, "data": master_data}

    log(' Done !')

    return assembled_data


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

def read_win_paths(file_paths, channel_table_path, fill_value=None, verbose=True):
    """
    Read WIN filepaths(s) into an ObsPy Stream.

    :param file_paths (str or list): Path or list of paths to WIN file(s)
    :param channel_table_path (str): Path to channel table file
    :param fill_value (float, optional): Fill value used in stream.merge() to handle gaps
    :param verbose (bool): If `False`, all print statements will be blocked. Default is `True`.
    :return: :class:`~obspy.core.stream.Stream`
    """

    # Use print for logging if verbose; otherwise use a no-op function
    log = print if verbose else lambda *args, **kwargs: None
    run_clock = time.time()

    # Load channel table
    log(f'Loading channel table...')
    channel_table = read_channel_table(channel_table_path)

    # Ensure file paths is a list
    if isinstance(file_paths, str):
        file_paths = [file_paths]
    # Check emptiness
    if len(file_paths) == 0:
        raise FileNotFoundError('Input list of files is empty.')
    log(f'Reading {len(file_paths)} file(s) to Stream...')

    # Initialize Stream to load data into
    stream = Stream()

    # Load each file in the list one by one
    for file_path in file_paths:
        assembled_data = _read_single_win_file(file_path, verbose)
        for sensor_id, sensor_data in assembled_data.items():
            row = channel_table[channel_table["location"] == sensor_id].iloc[0]
            network, station = row["station"].split(".")
            stats = Stats({
                "network": network,
                "station": station,
                "channel": row["channel"],  # adjust if column name differs
                "location": sensor_id,
                "sampling_rate": sensor_data["sampling_rate"],
                "starttime": sensor_data["starttime"],
                "npts": len(sensor_data["data"]),
            })
            trace = Trace(data=sensor_data["data"] * row["amplitude_correction"], header=stats)
            stream.append(trace)

    # Could also merge only once here... Seems slightly faster
    log('-----All files successfully read. Merging all streams...')
    _safe_merge(stream, fill_value)
    log('-----Streams merged.')

    log('Returning Stream object. Time elapsed: %s seconds.\n' % (time.time() - run_clock))

    return stream


def read_win(starttime, endtime, file_directory, file_pattern, file_interval, fill_value, channel_table_path, verbose=True):
    """
    Read WIN format data from a directory or archive, and trim

    :param starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time of desired data
    :param endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): End time of desired data
    :param file_directory (str): Path to directory containing WIN data
    :param file_pattern (str): WIN data file pattern (must be a valid strftime input for UTCDateTime)
    :param file_interval (str): WIN data file interval ("minute", "hour", or "day")
    :param fill_value (float): Value passed to stream.merge() when combining data from separate WIN files
    :param channel_table_path (str): Path to channel table file
    :param verbose (bool): If `False`, all print statements will be blocked. Default is `True`.
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
        file_times = np.arange(first_file_time, end_file_time, file_interval_time)
    else:
        file_times = np.arange(first_file_time, end_file_time + file_interval_time, file_interval_time)

    # Get list of WIN file paths
    file_paths = []
    for file_time in file_times:
        pattern = os.path.join(file_directory, file_time.strftime(file_pattern))
        matching_files = glob(pattern)
        file_paths.extend(matching_files)

    # Read WIN and trim
    stream = read_win_paths(file_paths, channel_table_path, fill_value=fill_value, verbose=verbose)
    stream = stream.trim(starttime, endtime)

    return stream


if __name__ == '__main__':

    starttime = UTCDateTime(2026, 3, 5, 10)  # UTC
    endtime = UTCDateTime(2026, 3, 5, 20)  # UTC
    file_directory = "./kurokami123/*/"  # for glob
    file_pattern = "%y%m%d%H"  # UTCDateTime strftime input
    file_interval = "hour"  # "minute", "hour", or "day"
    fill_value = None  # for stream.merge
    channel_table_path = "./channels.tbl"
    verbose = True

    stream = read_win(starttime,
                      endtime,
                      file_directory,
                      file_pattern,
                      file_interval,
                      fill_value,
                      channel_table_path,
                      verbose)
import struct
import numpy as np
import warnings
from obspy import UTCDateTime
from .utils import _to_signed_array, _uint2int

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

        # Conduct a check to see if the length of unit bytes matches header info
        self.read_state = "healthy"

        # Sensor ID and sampling rate is known, but bytes length is mismatched
        if len(unit_bytes) < self.size_bytes:
            warnings.formatwarning = lambda msg, cat, fname, lineno, line=None: f"{cat.__name__}: {msg}\n"
            warnings.warn("\nMismatch in number of bytes between unit header and channel block (%d/%d bytes). Falling back to sequential reading." %
                          (len(unit_bytes), self.size_bytes), RuntimeWarning)
            self.read_state = "truncated"

        # Reading channel blocks until the end of the unit
        while self.cursor_position < self.size_bytes:

            # Grab channel bytes
            channel_bytes = unit_bytes[self.cursor_position:self.size_bytes]

            # Insufficient bytes to derive sensor ID and sampling rate
            if len(channel_bytes) <= 4:
                warnings.warn("Insufficient bytes to derive sensor ID and sampling rate. Breaking.", RuntimeWarning)
                self.read_state = "fatal"
                break

            # Read channel block
            self._read_channel_block(channel_bytes)

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

        # Format strings controls how bytes are unpacked. Correct letter will be selected based on datasize
        format_string = 'BBHBI'

        # If read state is healthy, proceed
        if self.read_state == "healthy":

            # First data sample, always 4 bytes, in absolute units
            first_sample_value = _uint2int(struct.unpack(">I", channel_bytes[4:8])[0], 4)

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
            diff_data = _to_signed_array(diff_udata, datasize)

            # Convert to absolute units
            data_segment = np.cumsum(np.insert(diff_data, 0, first_sample_value)).astype(float)

        # Otherwise, read sequentially where possible, then pad with np.nan's
        elif self.read_state == "truncated":

            warnings.warn("Starttime %s; Sensor %s -- Padding with np.nan's." %
                          (self.starttime.strftime("%Y-%m-%d %H:%M:%S"), sensor_id), RuntimeWarning)

            # If even the first absolute sample is incomplete, return all nan's
            if len(channel_bytes) < 8:
                data_segment = np.full((sampling_rate,), np.nan)

            else:
                # First data sample, always 4 bytes, in absolute units
                first_sample_value = _uint2int(struct.unpack(">I", channel_bytes[4:8])[0], 4)

                # Read remaining data carefully, one sample at a time where possible
                diff_udata_list = []
                in_block_cursor_position = 8

                if datasize_code == 0:  # 0.5-byte packed data
                    while len(diff_udata_list) < (sampling_rate - 1):
                        if in_block_cursor_position >= len(channel_bytes):
                            break
                        current_byte = channel_bytes[in_block_cursor_position]
                        first_4bit = current_byte >> 4
                        second_4bit = current_byte & 0xf
                        diff_udata_list.append(first_4bit)
                        if len(diff_udata_list) < (sampling_rate - 1):
                            diff_udata_list.append(second_4bit)
                        in_block_cursor_position += 1

                elif datasize_code == 3:  # 3-byte packed data
                    while len(diff_udata_list) < (sampling_rate - 1):
                        if (in_block_cursor_position + 3) > len(channel_bytes):
                            break
                        first_8bit, second_8bit, third_8bit = struct.unpack(">BBB", channel_bytes[in_block_cursor_position:in_block_cursor_position + 3])
                        current_value = (first_8bit << 16) + (second_8bit << 8) + third_8bit
                        diff_udata_list.append(current_value)
                        in_block_cursor_position += 3

                else:
                    while len(diff_udata_list) < (sampling_rate - 1):
                        if (in_block_cursor_position + datasize) > len(channel_bytes):
                            break
                        current_value = struct.unpack(">" + format_string[datasize_code], channel_bytes[in_block_cursor_position:in_block_cursor_position + datasize])[0]
                        diff_udata_list.append(current_value)
                        in_block_cursor_position += datasize

                diff_udata = np.array(diff_udata_list)

                # Convert from unsigned to signed integers, still differential (i.e. difference to previous value rather than absolute value)
                diff_data = _to_signed_array(diff_udata, datasize)

                # Convert recovered portion to absolute units
                data_segment = np.cumsum(np.insert(diff_data, 0, first_sample_value)).astype(float)

                # Pad remaining samples with np.nan's
                if len(data_segment) < sampling_rate:
                    data_segment = np.concatenate((data_segment, np.full((sampling_rate - len(data_segment),), np.nan)))

        # Append data segment
        # Note that the data is still in integers, and not in physical units yet
        # Data will be converted to physical units all at once, after all 1s units are loaded
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


import warnings
import numpy as np

def _to_signed_array(uarray, byte_size):
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

def _uint2int(unsigned, byte_size):
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

def _get_opposite(number, byte_size):
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


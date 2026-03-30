import os
import time
import warnings
import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path
from obspy.core.trace import Stats
from obspy import Stream, Trace, UTCDateTime
from .onesecunit import OneSecUnit
from .utils import _safe_merge

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
# Import dependencies
from obspy import UTCDateTime
from read_win import read_win

# Read WIN data from directory
starttime = UTCDateTime(2026, 3, 5, 5)  # UTC
endtime = UTCDateTime(2026, 3, 5, 8)  # UTC
file_directory = "./read_win/sample_data/win_data"  # for glob
file_pattern = "%y%m%d%H"  # UTCDateTime strftime input
file_interval = "hour"  # "minute", "hour", or "day"
fill_value = None  # for stream.merge
channel_table_path = "./read_win/sample_data/channels.tbl"
verbose = True

stream = read_win(starttime,
                  endtime,
                  file_directory,
                  file_pattern,
                  file_interval,
                  fill_value,
                  channel_table_path,
                  verbose)

# Read bad WIN data to test warning messages
starttime = UTCDateTime(2026, 4, 2, 2)  # UTC
endtime = UTCDateTime(2026, 4, 2, 3)  # UTC

stream_bad = read_win(starttime,
                      endtime,
                      file_directory,
                      file_pattern,
                      file_interval,
                      fill_value,
                      channel_table_path,
                      verbose)
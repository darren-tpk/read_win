from obspy import UTCDateTime
from read_win import read_win

starttime = UTCDateTime(2026, 3, 5, 5)  # UTC
endtime = UTCDateTime(2026, 3, 5, 8)  # UTC
file_directory = "./read_win/sample_data"  # for glob
file_pattern = "%y%m%d%H"  # UTCDateTime strftime input
file_interval = "hour"  # "minute", "hour", or "day"
fill_value = None  # for stream.merge
channel_table_path = "./read_win/sample_data/channels.tbl"
verbose = False

stream = read_win(starttime,
                  endtime,
                  file_directory,
                  file_pattern,
                  file_interval,
                  fill_value,
                  channel_table_path,
                  verbose)
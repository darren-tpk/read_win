import os
import glob
import argparse
from pathlib import Path
from read_win import read_win_paths

# This script converts a list of WIN files into MSEED (one MSEED per channel).
# Usage:
# python convert_win2mseed.py /path/to/WIN_folder/ /path/to/mseed_folder/ /path/to/channel_table (--utc-offset 9 --fill-value 0 -v )

def main():

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Convert WIN files in a directory to miniSEED files.",
        allow_abbrev=False,
    )
    parser.add_argument("input_dir", help="Directory containing WIN files")
    parser.add_argument("output_dir", help="Directory for output miniSEED files")
    parser.add_argument("channel_table_path", help="Path to channel table file")
    parser.add_argument("--fill-value", default=0, type=float, dest="fill_value", help="Fill value for gaps in merged data")
    parser.add_argument("--utc-offset", default=0, type=float, dest="utc_offset", help="Deviation of WIN-file times from UTC, in hours (e.g., for JST, it is +9 from UTC)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print progress messages")
    input_args = parser.parse_args()

    # Check all input paths exist
    if not os.path.isdir(input_args.input_dir):
        raise NotADirectoryError(f"Input directory '{input_args.input_dir}' does not exist.")

    if not os.path.isdir(input_args.output_dir):
        raise NotADirectoryError(f"Output directory '{input_args.output_dir}' does not exist.")

    if not os.path.isfile(input_args.channel_table_path):
        raise FileNotFoundError(f"Channel table path '{input_args.channel_table_path}' does not exist.")

    list_files = glob.glob(os.path.join(input_args.input_dir, "*"))
    list_files.sort()

    # Make sure there are files to convert
    if len(list_files) == 0:
        raise FileNotFoundError("No WIN files found.")

    print("------------------------------------------------------------------")
    print("Beginning conversion process...")
    print("------------------------------------------------------------------")
    print(f"Found {len(list_files)} WIN file(s) to be converted.\n")

    # Loop through all files and convert them to a Stream one by one, then save each Trace in the Stream to MSEED
    for file_path in list_files:
        if input_args.verbose:
            print(f"Converting {file_path}")

        stream = read_win_paths(
            file_paths=[file_path],
            channel_table_path=input_args.channel_table_path,
            utc_offset=input_args.utc_offset,
            fill_value=input_args.fill_value,
            verbose=input_args.verbose,
        )

        # Loop through all Traces in the Stream and save to MSEED
        for trace in stream:
            trace.stats.location = str(trace.stats.location)[-2:] # MSEED only accepts location names with 2 characters max

            mseed_name = (
                f"{trace.stats.network}."
                f"{trace.stats.station}."
                f"{trace.stats.location}."
                f"{trace.stats.channel}."
                f"{trace.stats.starttime.year}."
                f"{trace.stats.starttime.julday:03d}."
                f"{trace.stats.starttime.hour:02d}"
            )

            save_filepath = Path(input_args.output_dir) / mseed_name
            trace.write(str(save_filepath), format="MSEED")

    print("------------------------------------------------------------------")
    print("...finished conversion to MSEED process.")
    print("------------------------------------------------------------------")


if __name__ == "__main__":
    main()

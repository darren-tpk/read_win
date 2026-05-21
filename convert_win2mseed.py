import os
import subprocess
import glob
import re
import json
import argparse
import obspy
import numpy as np
import matplotlib.pyplot as plt
import warnings
from pathlib import Path

import read_win

# This script converts WIN data files into MSEED files
# Converts all WIN files inside input folder
# How to use:
# source correct python environment
# python convert_win2mseed.py /path/to/WIN_folder/ /path/to/mseed_folder/ SA KUR (--fill-value 0) (-v) 


def main():

    # Set up command-line interface
    parser = argparse.ArgumentParser(description='Convert DATA-CUBE files to '
                                                 'miniSEED files while trimming, '
                                                 'adding metadata, and renaming. '
                                                 'Optionally extract coordinates '
                                                 'from digitizer GPS.',
                                     allow_abbrev=False)
    parser.add_argument('input_dir',
                        help='input directory containing WIN files ')
    parser.add_argument('output_dir',
                        help='directory for output miniSEED files')
    parser.add_argument('channel_table',
                        help='path to channel table')
    parser.add_argument('network',
                        help='desired SEED network code (2 characters, A-Z)')
    parser.add_argument('station',
                        help='desired SEED station code (3-4 characters, A-Z & '
                             '0-9)')
    parser.add_argument('--fill-value', default=0, type=float,
                        dest='fill_value',
                        help='fill_value for gaps in data')
    input_args = parser.parse_args()

    # Check if input directory is valid
    if not os.path.exists(input_args.input_dir):
        raise NotADirectoryError(f'Input directory \'{input_args.input_dir}\' doesn\'t '
                                     'exist.')

    # Check if output directory is valid
    if not os.path.exists(input_args.output_dir):
        raise NotADirectoryError(f'Output directory \'{input_args.output_dir}\' '
                                 'doesn\'t exist.')

    # Check network code format
    input_args.network = input_args.network.upper()
    if not re.fullmatch('[A-Z]{2}', input_args.network):
        raise ValueError(f'Network code \'{input_args.network}\' is not valid.')

    # Check station code format
    input_args.station = input_args.station.upper()
    if not re.fullmatch('[A-Z0-9]{3,4}', input_args.station):
        raise ValueError(f'Station code \'{input_args.station}\' is not valid.')

    # Find directory containing this script
    script_dir = os.path.dirname(__file__)

    print('------------------------------------------------------------------')
    print('Beginning conversion process...')
    print('------------------------------------------------------------------')

    # Print requested metadata
    print(f' Network code: {input_args.network}')
    print(f' Station code: {input_args.station}')


    # List all files in directory
    list_files = glob.glob( os.path.join(input_args.input_dir, '*') )
    list_files.sort()


    if len(list_files) == 0:
        raise FileNotFoundError('No WIN files found.')
    else:
        print(f'Found {len(list_files)} win files to be converted.\n')

    # Faster to do file one by one
    for file in list_files:

        # Read WIN data into Stream object
        stream = read_win.read_win(file, input_args.channel_table, network=input_args.network, station=input_args.station, fill_value=input_args.fill_value)

        # Rename location to match SEED convention
        for tr in stream:
            tr.stats.location = tr.stats.location[-2:]

        # Save each Trace object as MSEED
            mseed_name = f'{input_args.network}.{input_args.station}.{tr.stats.location}.{tr.stats.channel}.{tr.stats.starttime.year}.{tr.stats.starttime.julday}.{tr.stats.starttime.hour}'
            save_filepath = Path(input_args.output_dir) / mseed_name
            tr.write(save_filepath, format='MSEED')


    print('------------------------------------------------------------------')
    print('...finished conversion to MSEED process.')
    print('------------------------------------------------------------------')


if __name__ == '__main__':
    main()











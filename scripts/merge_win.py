import argparse

# This script merges WIN files into a single one.
# Usage:
# python merge_win.py /path/to/input_files* /path/to/output (-v)

def main():

    # Parse arguments
    parser = argparse.ArgumentParser(
                description="Merge WIN files."
            )
    parser.add_argument("input_files", nargs='+', help="List of WIN files to be merged [wildcards (*) supported]")
    parser.add_argument("output_file", help="Output WIN file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print progress messages")
    input_args = parser.parse_args()

    # Use print for logging if verbose; otherwise use a no-op function
    log = print if input_args.verbose else lambda *args, **kwargs: None

    log(f'Merging {len(input_args.input_files)} WIN files...')

    # Initialize variable for bytes to be read
    list_bytes = []

    # Cycle through all files and read all bytes
    for file in input_args.input_files:  
        with open(file,'rb') as f:
            list_bytes.append(f.read())

    # Merge all bytes into single variable
    merged_bytes = b''.join(list_bytes)

    # Save merged WIN file
    with open(input_args.output_file, "wb") as output_file:
        output_file.write(merged_bytes)

    log(f'Files successfully merged to {input_args.output_file}')

if __name__ == '__main__':
    main()

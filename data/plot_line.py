import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def plot_line(rrcs_list, window, outpath):
    """
    Plot RRCS values over time with a moving average.

    Parameters:
        rrcs_list (list of tuples): List of (frame, RRCS) values.
        window (int): Window size for moving average.
    """
    frames = [int(x[0]) for x in rrcs_list]
    rrcs_values = [x[1] for x in rrcs_list]

    if len(rrcs_values) < window:
        raise ValueError("Window size is larger than the available data points.")

    # Compute moving average
    moving_avg = np.convolve(rrcs_values, np.ones(window) / window, mode='valid')
    time_avg = frames[window - 1:]

    # Plot the RRCS values
    plt.figure(figsize=(8, 3))
    plt.plot(frames, rrcs_values, color='black', alpha=0.25, linewidth=0.75, label='RRCS')
    plt.plot(time_avg, moving_avg, color='red', linewidth=1, alpha=1, label=f'{window}-point Moving Average')

    plt.xlabel("Frame")
    plt.ylabel("RRCS")
    plt.legend()
    #plt.title("RRCS Time Evolution")
    plt.savefig(outpath, bbox_inches='tight', dpi=300)


def parse_file(infile):
    """
    Parse the RRCS output file.

    Parameters:
        infile (str): Path to the input file.

    Returns:
        dict: Dictionary with residue pairs as keys and lists of (frame, RRCS) tuples as values.
    """
    rrcs_dict = defaultdict(list)

    with open(infile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Frame"):
                continue
            frame, res1, res2, rrcs = line.split()
            rrcs_dict[f"{res1}-{res2}"].append((frame, round(float(rrcs), 3)))

    return rrcs_dict


def parse_arguments():
    """
    Parse command-line arguments using argparse.

    Returns:
        Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process RRCS output file and plot specified residue pair.")

    parser.add_argument("infile", type=str, help="Input file containing RRCS data.")
    parser.add_argument("outfile", type=str, help="Output file name (not used in this version).")
    parser.add_argument("pair", type=str, help="Residue pair to plot (e.g., 'A:ARG12-B:GLU34').")
    parser.add_argument("-w", "--window", type=int, default=25, help="Moving average window size (default: 25).")

    return parser.parse_args()


def main():
    """
    Main function to parse the RRCS file and generate a plot.
    """
    args = parse_arguments()
    rrcs_dict = parse_file(args.infile)

    if args.pair not in rrcs_dict:
        raise ValueError(f"Residue pair {args.pair} not found in the data.")

    plot_line(rrcs_dict[args.pair], args.window, args.outfile)


if __name__ == "__main__":
    main()

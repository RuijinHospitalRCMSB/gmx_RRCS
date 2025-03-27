import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def plot_bar_and_scatter(raw_data, residues, outpath):
    """
    Plot a bar chart for mean values and a scatter plot for individual data points.

    Parameters:
        raw_data (dict): Dictionary with keys as residue IDs and values as lists of RRCS values.
        residues (str): Comma-separated residue IDs to plot.
        outpath (str): Path to save the output plot.
    """
    # Clean and parse residue list
    residues = residues.replace("'", "").replace('"', "").strip()
    residue_list = [r.strip() for r in residues.split(",")]

    # Filter data for specified residues
    data = {key: raw_data[key] for key in residue_list if key in raw_data}

    if not data:
        print(f"Warning: None of the specified residues ({residues}) exist in the dataset.")
        return

    categories = list(data.keys())
    values = [data[key] for key in categories]

    # Compute mean RRCS values
    means = [np.mean(v) if v else 0 for v in values]

    # Create figure
    plt.figure(figsize=(12, 6))

    # Scatter plot with jitter
    all_x, all_y = [], []
    for i, (cat, v) in enumerate(zip(categories, values)):
        x_jitter = np.random.normal(loc=i, scale=0.05, size=len(v))  # Add slight jitter to avoid overlapping
        all_x.extend(x_jitter)
        all_y.extend(v)

    # Bar chart for mean values
    plt.bar(categories, means, color='lightblue', alpha=0.7, label='Mean RRCS')

    # Scatter points
    plt.scatter(all_x, all_y, alpha=0.15, color='black', s=5, label='Data Points')

    # Beautify plot
    plt.xlabel("Residue ID")
    plt.ylabel("RRCS")
    plt.grid(axis='y', linestyle='--', alpha=0.6)

    # Save figure
    plt.savefig(outpath, bbox_inches='tight', dpi=300)
    print(f"Plot saved to {outpath}")


def parse_file(infile):
    """
    Parse the RRCS output file.

    Parameters:
        infile (str): Path to the input file.

    Returns:
        dict: Dictionary with residue IDs as keys and lists of RRCS values.
    """
    rrcs_dict = defaultdict(list)
    pattern = re.compile(r"(\w{1,8}):(\d{1,5})(\w{3})")  # Regex pattern for extracting residue ID

    with open(infile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Frame"):
                continue
            parts = line.split()
            if len(parts) != 4:
                print(f"Skipping malformed line: {line}")
                continue

            frame, res1, res2, rrcs = parts

            # Extract residue ID using regex
            match = pattern.match(res1)
            if not match:
                print(f"Skipping invalid residue format: {res1}")
                continue

            res1_id = match.group(2)  # Extract numerical residue ID
            try:
                rrcs_value = round(float(rrcs), 3)
                rrcs_dict[res1_id].append(rrcs_value)
            except ValueError:
                print(f"Skipping invalid RRCS value: {rrcs}")
                continue

    return rrcs_dict


def parse_arguments():
    """
    Parse command-line arguments using argparse.

    Returns:
        Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process RRCS output file and plot specified residues.")

    parser.add_argument("-i", "--infile", type=str, required=True, help="Input file containing RRCS data.")
    parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file for the generated plot.")
    parser.add_argument("-r", "--residues", type=str, required=True,
                        help="Comma-separated residue IDs to plot (e.g., '25,26,27').")

    return parser.parse_args()


def main():
    """
    Main function to parse the RRCS file and generate a plot.
    """
    args = parse_arguments()
    rrcs_dict = parse_file(args.infile)
    plot_bar_and_scatter(rrcs_dict, args.residues, args.outfile)


if __name__ == "__main__":
    main()

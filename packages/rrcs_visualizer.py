import sys
import argparse

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

try:
    # Attempt to import modules from a subdirectory named 'packages' within the current directory.
    from packages.utilities import *
except:
    try:
        # If the initial import fails (likely due to differing working directories), 
        # attempt to import using a relative path.
        from .utilities import *
    except:
        # If importing with a relative path also fails, 
        # attempt to import directly from the current directory.
        from utilities import *

def parse_command():
    """config arguments from command line"""
    parser = argparse.ArgumentParser(description="Process configuration parameters.")
    
    # define required parameter
    parser.add_argument('--top_file', type=str, required=True, 
                                help="Topology file path (required)")
    parser.add_argument('--traj_file', type=str, required=True, 
                                help="Trajectory file path (required)")
    
    # define optional parameter
    parser.add_argument('--res_file', type=str, default='',
                                help="Path to the file containing residue pair indices.")

    if len(sys.argv) == 1:
        log_error("InvalidParameter", "No arguments provided. Displaying help:\n")
        self.parser.print_help()
        sys.exit(1)
    args = self.parser.parse_args()
    return vars(args)
    


def plot_line():
    # Generate some sample data
    np.random.seed(0)
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x) + np.random.normal(0, 0.1, len(x))
    # y2 = np.cos(x) + np.random.normal(0, 0.1, len(x))

    # Create a seaborn style plot
    sns.set(style="whitegrid")

    # Initialize the matplotlib figure
    plt.figure(figsize=(12, 6))

    # Plot the data
    plt.plot(x, y1, label='Sine', color='b', linestyle='-', linewidth=1, marker='o')
    # plt.plot(x, y2, label='Cosine', color='r', linestyle='--', linewidth=2, marker='x')

    # Add a title and labels
    plt.title('Sine and Cosine Waves with Noise', fontsize=20, fontweight='bold')
    plt.xlabel('X-axis', fontsize=14)
    plt.ylabel('Y-axis', fontsize=14)

    # Customize the ticks on the axes
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Add a legend
    plt.legend(loc='upper right', fontsize=12)

    # Remove the top and right spines for a cleaner look
    sns.despine()

    # Show the plot
    plt.show()


def main():
    plot_line()


if __name__ == "__main__":
    main()

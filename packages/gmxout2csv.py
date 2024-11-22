# -*- coding: utf-8 -*-

"""
This script converts text files output by GMX_RRCS to CSV format for easier subsequent processing,
and checks that all rows have the same number of elements.
"""

import sys
import argparse
from collections import Counter

def transform_text_to_csv(txt_file):
    """
    Transforms a text file into a CSV format.
    
    Args:
        txt_file (str): The path to the input text file.
        
    Returns:
        list: A list of strings where each string represents a line in CSV format.
    """
    outlines = []
    n_elements = []
    line_count = 0
    with open(txt_file, 'r') as file:
        for line in file:
            line_count += 1
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            else:
                elements = line.split()
                outlines.append(','.join(elements))
                n_elements.append((len(elements), line_count))
    return outlines, n_elements

def check_row_consistency(n_elements):
    """
    Checks if all rows have the same number of elements.
    
    Args:
        n_elements (list): A list containing the number of elements in each row.
        
    Returns:
        bool: True if all rows have the same number of elements, False otherwise.
    """
    n_elements_per_line = [n for n, _ in n_elements]
    stand = sorted(Counter(n_elements_per_line).items(), key=lambda x: x[1], reverse=True)[0][0]
    if len(set(n_elements_per_line)) > 1:
        for n_elements, line_number in n_elements:
            if n_elements != stand:
                print(f"Line {line_number} has {n_elements} elements. check the file")
        return False
    return True

def write_to_csv(outlines, csv_file):
    """
    Writes the transformed data to a CSV file.
    
    Args:
        outlines (list): A list of strings in CSV format.
        csv_file (str): The path to the output CSV file.
    """
    with open(csv_file, 'w') as file:
        file.writelines([f"{line}\n" for line in outlines])

def parase_arguments():
    """
    Parse the command line arguments.
    
    This function uses the argparse module to define and parse the command line arguments needed by the program.
    It includes the input file path and the output file path, both of which are optional.
    
    Returns:
        The parsed command line arguments as a Namespace object.
    """
    parser = argparse.ArgumentParser(description="Transform GMX_RRCS output file to CSV format.")
    parser.add_argument('-i', '--input', type=str, help='Path to the input text file.')
    parser.add_argument('-o', '--output', type=str, help='Path to the output csv file.')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    return args


def main():
    # if len(sys.argv) != 3:
    #     print("Usage: python script.py <input_txt_file> <output_csv_file>")
    #     sys.exit(1)
    
    args = parase_arguments()
    
    txt_file = args.input
    csv_file = args.output
    
    outlines, n_elements = transform_text_to_csv(txt_file)
    if check_row_consistency(n_elements):
        write_to_csv(outlines, csv_file)
    else:
        print("Data not written to CSV due to inconsistency in row lengths.")


if __name__ == '__main__':
    main()

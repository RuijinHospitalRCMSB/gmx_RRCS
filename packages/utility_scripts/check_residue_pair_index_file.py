"""
Script Functionality Description

This Python script validates text files that describe 
a custom system of residue index pairs. Its key features are:

1. Format Validation:
   - Checks every line against predefined rules.
   - Ensures absence of illegal characters except digits,
     spaces, hyphens '-', and dollar signs '$'.
   - Verifies correct counts of '$' symbols per line.
   - Parses numerical ranges within lines.

2. Error Detection & Reporting:
   - Captures deviations from format rules.
   - Lists detailed information about format violations.
   - Includes line numbers and error descriptions.

3. Data Parsing:
   - Converts numerical ranges into lists of numbers.
   - Transforms lines into structured formats.

4. User-Friendly Execution:
   - Accepts file path as an argument.
   - Enables direct command-line execution.

By integrating these features, the script ensures accuracy 
and consistency in residue index pair files. It prevents 
data analysis issues due to formatting errors, making it 
suitable for scientific research, bioinformatics, and 
fields requiring meticulous data quality control.
"""


import re
import sys

def parse_line(line):
    # Remove comments and strip whitespace
    line = re.sub(r";.*", "", line).strip()
    
    # Check for illegal characters
    if re.search(r"[^0-9 \-$]", line):
        return False, "Line contains illegal characters."
    
    # Check if the line contains the '$' separator
    if "$" in line:
        parts = line.split("$")
        if len(parts) != 2:
            return False, "Line contains more than one '$' symbol."
        
        left, right = map(lambda x: x.strip().split(), parts)
        left = [parse_range(item) for item in left]
        right = [parse_range(item) for item in right]
        return True, (left, right)
    else:
        items = line.split()
        items = [parse_range(item) for item in items]
        return True, items

def parse_range(item):
    if "-" in item:
        start, end = map(int, item.split("-"))
        return list(range(start, end + 1))
    else:
        return [int(item)]

def validate_file(file_path):
    errors = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            valid, result = parse_line(line)
            if not valid:
                errors.append(f"Error at line {i+1}: {result}")
    
    if errors:
        print("Validation errors:")
        for error in errors:
            print(error)
    else:
        print("File is correctly formatted.")


if __name__ == "__main__":
    infile = sys.argv[1] # your input file path
    validate_file(infile)


# Example
# python check_residue_pair_index_file.py format_test\SpecifyResiduePairIndex2.txt
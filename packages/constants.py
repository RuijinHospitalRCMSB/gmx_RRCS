# -*- coding: utf-8 -*-

"""
Global constants used across different modules.

This module contains various constants used throughout the application,
including mappings, thresholds, and lists of specific types.
"""


# Mapping dictionary from three-letter to one-letter amino acids
THREE_TO_ONE_LETTER = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X'
}

# Protein main chain atom types
MAIN_CHAINS_ATOMS = ['N', 'CA', 'C', 'O']

# Define distance threshold as a constant for better readability
ATOM_DISTANCE_THRESHOLD = 4.14

# Maximum index difference for neighbors
MAX_INDEX_DIFFERENCE_FOR_NEIGHBORS = 5

# The minimum amount of space between columns in the output file.
OUTPUT_OFFSET = 2

# The minimum amount of frames to be considered for the calculation
STEP_INDEX_LOWER_LIMIT = 5

# Warning threshold for excessive trajectory computation frames
STEP_INDEX_WARN_LIMIT = 10000

# the expression used to exclude all hydrogen atom names
prefixs = ['', '1', '2', '3']
AND_NOT_H_ATOMS = ""
for prefix in prefixs:
    AND_NOT_H_ATOMS += f"and not name {prefix}H* "


# Header description of the output file
OUTPUT_HEADER = """
# This document provides a detailed record of the calculation results of the 
# Residue-Residue Contact Score (RRCS) from GROMACS trajectories. RRCS is a 
# quantitative metric used to describe the interactions between residues within 
# a protein, specifically reflecting the degree of proximity between residues 
# in three-dimensional space. This metric is crucial for understanding the 
# dynamic changes in protein structures and their functional realization.

# Columns:
# - Frame: Indicates the frame index in the GROMACS simulation trajectory.
# - Residue1: Lists the identifier of the first residue involved in the RRCS calculation.
# - Residue2: Lists the identifier of the second residue involved in the RRCS calculation.
# - RRCS: Displays the RRCS calculation score for the corresponding residue pair.

# Note: For a comprehensive understanding of RRCS theory and its critical applications 
# in protein structure analysis, please refer to our research published in 
# eLife, "Common activation mechanism of class A GPCRs." This study explains 
# the principles of RRCS and demonstrates its use in revealing the activation 
# mechanisms of class A GPCRs.



"""


# Header description of the output file
OUTPUT_FILTER_HEADER = """
# This document provides a detailed record of the calculation results of the 
# Residue-Residue Contact Score (RRCS) from GROMACS trajectories. RRCS is a 
# quantitative metric used to describe the interactions between residues within 
# a protein, specifically reflecting the degree of proximity between residues 
# in three-dimensional space. This metric is crucial for understanding the 
# dynamic changes in protein structures and their functional realization.

# Columns:
# - Frame: Indicates the frame index in the GROMACS simulation trajectory.
# - Residue1: Lists the identifier of the first residue involved in the RRCS calculation.
# - Residue2: Lists the identifier of the second residue involved in the RRCS calculation.
# - RRCS: Displays the RRCS calculation score for the corresponding residue pair.

# Note: Here are the refined results, where only those with higher scores are retained.

# Note: For a comprehensive understanding of RRCS theory and its critical applications 
# in protein structure analysis, please refer to our research published in 
# eLife, "Common activation mechanism of class A GPCRs." This study explains 
# the principles of RRCS and demonstrates its use in revealing the activation 
# mechanisms of class A GPCRs.



"""
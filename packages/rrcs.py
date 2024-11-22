# -*- coding: utf-8 -*-

import os
import sys
import itertools
import logging
import argparse

import timeit
import numpy as np
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt

from numba import jit
from collections import defaultdict
from typing import Dict, Union
from termcolor import colored 
from MDAnalysis import Universe
from concurrent.futures import wait, ALL_COMPLETED
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    # Attempt to import modules from a subdirectory named 'packages' within the current directory.
    from packages.utilities import *
    from packages.constants import *
except:
    try:
        # If the initial import fails (likely due to differing working directories), 
        # attempt to import using a relative path.
        from .utilities import *
        from .constants import *
    except:
        # If importing with a relative path also fails, 
        # attempt to import directly from the current directory.
        from utilities import *
        from constants import *
finally:
    from utilities import *
    from constants import *


# configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class ConfigParser:
    def __init__(self):
        """config arguments from command line"""
        self.parser = argparse.ArgumentParser(description="Process configuration parameters.")
        
        # define required parameter
        self.parser.add_argument('--pdb_file', type=str, required=True, 
                                 help="PDB file path (required)")
        
        # define optional parameter
        self.parser.add_argument('--res_file', type=str, default='',
                                 help="Path to the file containing residue pair indices.")
        self.parser.add_argument('--output_file', type=str, default='RRCS_output.txt',
                                 help="Output file name, default is 'RRCS_output.txt'")
        self.parser.add_argument('--filter_threshold', type=float, default=3.0,
                                 help='Choose whether to output the high-scoring results that have been filtered, default is 3.0.')
        self.parser.add_argument('--num_processes', type=int, default=None,
                                 help="Number of processes for parallel execution."
                                  "If None, use all available CPU cores. Default is None.")
        self.parser.add_argument('--radius_min', type=float, default=3.23,
                                 help="Minimum distance threshold in Ångström, default is 3.23")
        self.parser.add_argument('--radius_max', type=float, default=4.63,
                                 help="Maximum distance threshold in Ångström, default is 4.63")


    def parse_arguments(self):
        """Parses command line arguments and returns them as a dictionary."""
        if len(sys.argv) == 1:
            log_error("InvalidParameter", "No arguments provided. Displaying help:\n")
            self.parser.print_help()
            sys.exit(1)
        args = self.parser.parse_args()
        return vars(args)
 

@timing_decorator
class ResidueCombinePairs:
    def __init__(self, basic_settings):
        """
        Initialize a ResidueCombinePairs object.
        
        Parameters:
            basic_settings (dict): A dictionary containing basic settings.
        """
        self.basic_settings = basic_settings
        self.res_pair_set = set()

    def read_file(self):
        """Read and parse the index file, populate res_pair_set."""
        res_file = self.basic_settings["res_file"]
        if os.path.exists(res_file):
            with open(res_file, 'r') as f:
                lines = f.readlines()
            lines = self.clear_lines(lines)
            if len(lines) == 0:
                log_warning(
                    "InputFileWarning", 
                    "The ResidueIndex file is empty. All residue pairs will be processed, which may take a considerable amount of time."
                    )
                self.res_pair_set = self.combine_all_residues()
            else:
                for line in lines:
                    self.res_pair_set.update(self.parse_line(line))
        else:
            log_warning(
                "InputFileWarning", 
                "The ResidueIndex file does not exist. All residue pairs will be processed, which may take a considerable amount of time."
                )
            self.res_pair_set = self.combine_all_residues()
            
    def clear_lines(self, lines):
        """Clear the lines in the ResidueIndex file."""
        new_lines = []
        for line in lines:
            line = line.strip()
            if line:  # Check if line is not empty
                new_lines.append(line)
        return new_lines

    def parse_line(self, line):
        """
        Uniformly parse a line, handling both sides of the '$' symbol similarly.
        
        Parameters:
            line (str): The line to be parsed.
            
        Returns:
            tuple: The result of the parsing.
        """
        line = line.strip().split(";")[0].strip() # Remove comments
        # If the line consists only of alphanumeric characters, it is processed as a special type of line
        if is_alnum_space(line):
            return self.parse_other_line(line)
        # If the line contains the '$' symbol, it is processed as a standard type of line
        elif '$' in line:
            return self.parse_stand_line(line)
        else:
            log_error(
                "InputFileError", 
                "The --res_file you entered contains illegal characters, please ensure that your delimiter is '$'."
                )
        
    def parse_stand_line(self, line):
        """Parse lines containing '$', generating all possible combinations of residue pairs."""
        res_pair_list = []
        parts = line.split("$")
        for part in parts:
            part = part.strip()  # Remove leading/trailing whitespaces
            if part:  # Check if part is not empty
                res_selc = self._parse_res_list(part, self.basic_settings['res_num'])
                res_pair_list.append(res_selc)
            else:
                raise ResidueIndexError
        assert len(res_pair_list) == 2, ResidueIndexError
        return set(itertools.product(*res_pair_list))
    
    def parse_other_line(self, line):
        """Parse lines not containing with '$', generating all possible combinations of residue pairs."""
        res_selc = self._parse_res_list(line.strip(), self.basic_settings['res_num'])
        return set(itertools.combinations(sorted(res_selc), 2))

    @staticmethod
    def _parse_res_list(res_list_str, res_num):
        """Parse a residue list from a string, supporting range notation."""
        res_selc = []
        for res in res_list_str.split():
            if '-' in res:
                start, end = map(int, res.split('-'))
                res_selc.extend(range(start, end + 1))
            else:
                res_selc.append(int(res))
        if 'all' in res_list_str:
            res_selc = list(range(1, res_num + 1))
        return set(res_selc)
    
    def combine_all_residues(self):
        """
        Generates all possible residue pairs from an MD trajectory.

        This function checks if the MD trajectory has a 'residues' attribute,
        extracts the residue IDs, sorts them if necessary, and then generates
        all unique pairs of consecutive residues. It handles edge cases such as
        missing keys or attributes in the input data and catches any unexpected
        exceptions to ensure robustness.

        Returns:
            generator: A generator that yields tuples representing each pair of
                    consecutive residues in the sorted list.

        Raises:
            KeyError: If 'protein' key is missing in the basic settings.
            AttributeError: If 'residues' attribute is missing in the MD trajectory.
            Exception: For any other unexpected errors during execution.
        """
        protein = self.basic_settings['protein']
        # Extract resid values from residues
        residues = [res.resid for res in protein.residues]
        return set(itertools.combinations(residues, r=2))

    def get_res_pairs(self):
        """Retrieve the parsed set of residue pairs."""
        pairs = tuple(sorted(self.res_pair_set))
        logging.info(f"Read {len(pairs)} residue pairs. Proceeding to calculate RRCS between these {len(pairs)} residue pairs.")
        return pairs


@timing_decorator
class UniverseInitializer:
    def __init__(self, basic: dict):
        self.basic = basic

    def check_radius(self):
        """
        Ensure the radius minimum and maximum values are within a valid range.
        
        This method adjusts the minimum and maximum radius values in the basic configuration to ensure they are not less than 0.
        If the minimum radius is less than 0, it is set to 0; if the maximum radius is less than the minimum, the maximum is then set to the minimum.
        This process guarantees the reasonableness of radius settings, avoiding invalid or unreasonable search ranges.
        
        Returns:
            No return value, but modifies the 'radius_min' and 'radius_max' values in the self.basic dictionary.
        """
        # Ensure the radius minimum is non-negative
        self.basic['radius_min'] = max(self.basic['radius_min'], 0)
        self.basic['radius_max'] = max(self.basic['radius_min'], self.basic['radius_max'])
        logging.info(f"Radius settings updated. Minimum radius: {self.basic['radius_min']}, Maximum radius: {self.basic['radius_max']}.")

    def check_file_exists(self):
        """
        Verifies the existence of required input files.
        
        This method iterates through a list of file paths, checking if each file exists at the corresponding path. 
        If a file is not found, an exception is raised.
        """
        pass


    def initialize_universe(self):
        """
        Initialize the Universe object.

        This method creates an MDAnalysis Universe object. It uses the topology file and trajectory file paths retrieved from the 
        `basic` dictionary to instantiate the Universe. Upon initialization, it stores the number of residues in the 
        universe back into the `basic` dictionary for future reference.

        Returns:
            None
        """
        # Instantiate the Universe object with topology and trajectory file paths
        self.basic['md_traj'] = Universe(self.basic['pdb_file'])
        # Store the number of residues in the universe into the basic dictionary
        self.basic['protein'] = self.basic['md_traj'].select_atoms('protein')
        self.basic['res_num'] = len(self.basic['protein'].residues)
        logging.info(f"Topology file contains {self.basic['res_num']} residues.")

        traj_time = self.basic['md_traj'].trajectory[-1].time
        logging.info(f'Trajectory file contains {traj_time} ps.')

        self.basic['n_traj_steps'] = self.basic['md_traj'].trajectory.n_frames
        logging.info(f"Trajectory file contains {self.basic['n_traj_steps']} frames.")

        parser = ResidueCombinePairs(self.basic)
        parser.read_file()
        self.basic["res_pairs"] = parser.get_res_pairs()

    def get_chain(self):
        """
        Retrieves the chain information from the trajectory.
        """
        chains = []
        for chain in self.basic["protein"].segments:
            chains.append((chain.segindex, chain.segid))
        self.basic['traj_chains'] = chains

    def run(self):
        """
        Conducts a series of checks to ensure all prerequisites for the program's execution are met.
        
        This method invokes multiple verification functions to respectively confirm the availability of the plotting library,
        the reasonableness of the set radius, the existence of the RRCs file, the presence of required files, 
        the initialization of the universe model, the appropriateness of the minimum time step, and the validity of the time interval settings.
        These checks guarantee that the program operates under the correct configurations and conditions, averting potential errors and exceptions.
        """
        
        # Check if the specified radius is reasonable
        self.check_radius()
        # Verify the existence of necessary files
        self.check_file_exists()
        # Initialize the universe model
        self.initialize_universe()
        # Get chains from trajectory
        self.get_chain()


@timing_decorator
class DataVisualizer:
    def __init__(self, basic_settings, rrcs_data):
        """
        Initializes the class with basic settings.

        This constructor takes a dictionary of basic settings, most importantly 
        the path to the output directory. If the output directory does not exist,
        it will be created.

        :param basic_settings: A dictionary containing basic settings, must include 
                            the path to the output directory.
               rrcs_data: A dictionary containing the RRCs data.
        """
        
        # Concatenate the output file path based on the output directory path
        self.output = basic_settings['output_file']
        # Save whether to generate plots from the settings
        self.rrcs_data = rrcs_data
        self.filter_threshold = basic_settings['filter_threshold']

    def write_output(self):
        """Writes the processed data to the output file."""
        # Initialize a list with the output file's header information.
        outlines = [OUTPUT_HEADER]
        filter_outlines = [OUTPUT_FILTER_HEADER]

        # Append the reformatted data lines.
        lines, filter_lines = self.reformat_data_lines(self.rrcs_data)

        outlines.extend(lines)
        with open(self.output, 'w') as f:
            f.writelines(outlines)
        logging.info(f"RRCS data is saved to {self.output} file.")

        # Construct the file path for the filtered data.
        filepath = f"{os.path.splitext(self.output)[0]}_filter_rrcs_greater_than_{self.filter_threshold}.txt"
        filter_outlines.extend(filter_lines)
        with open(filepath, 'w') as f:
            f.writelines(filter_outlines)
        logging.info(f"Filtered RRCS data is saved to the {filepath} file.")

    def reformat_data_lines(self, rrcs_data):
        """
        Reformat data lines to create a uniform representation for each pair of residues and their corresponding 
        Relative Rotamer Conformation Score (RRCS).
        
        This method starts by adding a header, then transforms the data within the rrcs_data dictionary into a 
        standardized list format, and finally prints these data in a neat table format through the pretty_print_table 
        method.
        
        Returns:
            Formatted data table, returned as a string.
        """
        outlines = [('Frame', 'Residue1', 'Residue2', 'RRCS')]
        filter_outlines = [('Frame', 'Residue1', 'Residue2', 'RRCS')]
        for frame, rrcs_list in sorted(rrcs_data.items(), key=lambda x: x[0]):
            # Iterate over the RRCS list, each entry is a tuple containing residue identifiers and an RRCS score
            for res1, res2, rrcs_score in rrcs_list:
                outlines.append((frame, res1, res2, rrcs_score))
                if rrcs_score > self.filter_threshold:
                    filter_outlines.append((frame, res1, res2, rrcs_score))
        outlines = self.pretty_print_table(outlines)
        filter_outlines = self.pretty_print_table(filter_outlines)

        return outlines, filter_outlines

    def pretty_print_table(self, rows):
        """
        Prints a table in a pretty format with specified column widths.
        
        :param rows: A list of lists where each inner list represents a row in the table.
        :param output_offset: An integer representing the offset to add to each column's width.
        :return: A string representation of the formatted table.
        """
        # Calculate column widths
        col_widths = [max(len(str(item))+OUTPUT_OFFSET for item in col) for col in zip(*rows)]
        # Builds the formatted table string based on the calculated column widths.
        table = [
            ''.join(str(item).ljust(width) for item, width in zip(row, col_widths))
            for row in rows
        ]
        return '\n'.join(table)

    def run(self):
        """
        Executes the analysis and outputs the result.

        If the plotting option is configured, it will draw a scatter plot first. 
        Afterwards, regardless of whether plotting is done or not, it will write to the output file.
        """
        # Write to the output file
        self.write_output()


class RRCSAnalyzer:
    def __init__(self):
        # self.basic_settings = basic_settings
        pass

    def get_residue_name(self, protein, index_i, chain_ix):
        """
        Get residue names from the universe.
        """
        atoms = list(protein.select_atoms(f"resid {index_i} and segindex {chain_ix} {AND_NOT_H_ATOMS}").ids - 1)
        if atoms:
            res_name = THREE_TO_ONE_LETTER.get(protein.atoms[atoms[0]].resname, 'X')
        else:
            res_name = 'X'
        return res_name


    def are_residues_adjacent(self, index_i, index_j):
        """
        Check if two residues are neighbors within a specified distance.
        
        Parameters:
        - index_i: Integer representing the index of the first residue.
        - index_j: Integer representing the index of the second residue.
        - distance_threshold: The maximum distance (in Angstroms) considered as 'neighbor'. Default is 5A.
        
        Returns:
        - Boolean: True if the residues are neighbors within the specified distance, False otherwise.
        """
        return abs(index_i - index_j) < MAX_INDEX_DIFFERENCE_FOR_NEIGHBORS


    @staticmethod
    def adjest_atom_coordinates(is_neighbor, residue, frame_step):
        """
        Adjust the coordinates of the atoms in the residue.
        
        Parameters:
        is_neighbor (bool): Whether the residues are neighbors.
        residue (list): List of tuples, each containing:
                    - atom_name (str): Name of the atom.
                    - atom_ids (list): Atom ids list.
                    - atom_occupancy (float): Occupancy of the atom.
        
        Returns:
        numpy.ndarray: Adjusted coordinates of the atoms.
        """
        try:
            positions = frame_step.positions
            adjest_coord = [
                positions[atom_id] * atom_occupancy
                for atom_name, atom_id, atom_occupancy in residue
                if not (is_neighbor and atom_name in MAIN_CHAINS_ATOMS)
            ]
            return np.array(adjest_coord)
        except (TypeError, ValueError) as e:
            log_warning(e, "Error processing atom coordinates")
            return np.array([])


    # @jit(nopython=True)
    @staticmethod
    def prefilter_contacts(coord_i, coord_j):
        """
        Pre-filter the contacts to reduce the number of calculations.
        
        Parameters:
        coord_i, coord_j: Numpy arrays representing coordinates. Both are expected to have shape (n, 3),
                            where n is the number of coordinates.
        
        Returns:
        Boolean value indicating if there is at least one pair of coordinates closer than 4.14 units.
        """
        # Check if coord_i and coord_j are empty
        if (not coord_i.size) or (not coord_j.size):
            return False
        
        # Compute coordinate differences and filter pairs closer than DISTANCE_THRESHOLD
        # coord_i shape: (m, 3) -> (m, 1, 3)
        # coord_j shape: (n, 3) -> (1, n, 3)
        # diff shape: (m, n, 3)
        diff = np.abs(coord_i[:, np.newaxis, :] - coord_j[np.newaxis, :, :])

        return np.any(np.all(diff < ATOM_DISTANCE_THRESHOLD, axis=2))


    @staticmethod
    @jit(nopython=True)
    def get_distances(coord_i, coord_j):
        """
        Calculate the distances between two sets of coordinates using scipy's cdist function.
        
        Args:
            coord_i: A 2D numpy array representing the first set of coordinates. Each row is a coordinate.
            coord_j: A 2D numpy array representing the second set of coordinates. Each row is a coordinate.
        
        Returns:
            A 2D numpy array where element (i,j) represents the distance between coord_i[i] and coord_j[j].
        """
        diff = coord_i[:, np.newaxis, :] - coord_j[np.newaxis, :, :]
        return np.sqrt(np.sum(diff**2, axis=-1))


    @staticmethod
    @jit(nopython=True)
    def compute_rrcs_jit(distances_matrix, d_max_squared, d_min_squared):
        """
        JIT-compiled function to calculate RRCS more efficiently.
        
        Parameters:
        - distances_matrix (np.ndarray): Matrix of squared distances.
        - d_max_squared (float): Squared maximum distance for contact consideration.
        - d_min_squared (float): Squared minimum distance for full score.
        
        Returns:
        - total_score (float): The computed RRCS.
        """
        # Apply conditions using NumPy's where function
        scores = np.where(
            distances_matrix >= d_max_squared,
            # When distance is greater than or equal to max squared distance, score is 0
            0.0,
            np.where(
                distances_matrix <= d_min_squared,
                # When distance is less than or equal to min squared distance, score is 1
                1.0,
                # Calculate score between the two distances
                (d_max_squared - distances_matrix) / (d_max_squared - d_min_squared)
            )
        )

        # Sum all elements in the score matrix
        return np.sum(scores)


    def get_residue_info(self, protein, chains, residues):
        """
        Retrieve residue information from the universe.
        
        Parameters:
        residues (list): List of residue IDs.
        
        Returns:
        - A dictionary containing residue information.
        """
        pair_chain = defaultdict(dict)
        for _ix, _id in chains:
            _id = 'A' if _id == 'SYSTEM' else _id
            pair_residue = defaultdict(dict)
            for resid in residues:
                atom_ids = list(protein.select_atoms(f"resid {resid} and segindex {_ix} {AND_NOT_H_ATOMS}").ids - 1)
                res_name = ''
                pair_atom = []
                for atom_id in atom_ids:
                    atom = protein.atoms[atom_id]
                    atom_name = atom.name
                    # For each atom, attempt to retrieve its occupancy value
                    try:
                        atom_occu = atom.occupancy
                    except mda.exceptions.NoDataError as e:
                        # If a NoDataError exception occurs, this means the atom lacks occupancy information
                        # In such cases, we set the occupancy atom_occu to a default value of 1.0
                        atom_occu = 1.0
                    res_name = atom.resname
                    pair_atom.append((atom_name, atom_id, atom_occu))
                if pair_atom or res_name:
                    res_name = THREE_TO_ONE_LETTER.get(res_name, 'X')
                    pair_residue[f'{resid}{res_name}'] = pair_atom
            pair_chain[(_ix, _id)] = pair_residue
        return pair_chain


    def load_residues(self, res_pairs):
        """
        Processes a given list of tuples, extracting the first and second elements 
        from each tuple into two separate sets to eliminate duplicates. Finally, 
        it returns two sorted lists containing these unique elements.

        Args:
        - res_pairs: A list containing tuples. The first and second elements from 
                    each tuple will be extracted.

        Returns:
        - Two sorted lists containing all unique first and second elements from 
        the tuples, respectively.
        """

        # Input validation
        if not isinstance(res_pairs, tuple):
            logging.error("res_pairs must be a tuple")
        if not all(isinstance(pair, tuple) and len(pair) == 2 for pair in res_pairs):
            logging.error("Each element in res_pairs must be a tuple of length 2")
        
        member_first = set()
        member_second = set()
        
        for pair in res_pairs:
            member_first.add(pair[0])
            member_second.add(pair[1])
        member_first = sorted(member_first)
        member_second = sorted(member_second)

        return member_first, member_second


    def analyze_frame(self, frame_index, info_first, info_second, settings, md_traj):
        """
        Analyzes a specified frame of molecular structure to calculate the distances 
        and Relative Residual Contact Scores (RRCS) between specific residue pairs.
        
        Parameters:
        info_first : dict
            A dictionary containing residue information for the first model.
        info_second : dict
            A dictionary containing residue information for the second model.
        settings : dict
            A dictionary containing calculation settings such as residue pairs,
            minimum and maximum radii.
        md_traj : object
            A molecular dynamics trajectory object that holds structural information
            for all frames.
        
        Returns:
        frame_count : int
            The count of the current frame.
        frame_rrcs : list
            A list of all calculated RRCS values for the current frame.
        """
        # Retrieve step information for the specified frame
        frame_step = md_traj.trajectory[frame_index]
        protein = md_traj.select_atoms("protein")
        frame_count = frame_step.frame + 1
        # Initialize the list for RRCS values of the current frame
        frame_rrcs = []

        # Iterate over all chains and residues information in the first model
        for chain_ix, chain_id in info_first.keys():
            info_res_first = info_first[(chain_ix, chain_id)]
            info_res_second = info_second[(chain_ix, chain_id)]
            # Iterate over the residue pairs defined in settings
            for index_i, index_j in settings['res_pairs']:
                res_i = self.get_residue_name(protein, index_i, chain_ix)
                res_j = self.get_residue_name(protein, index_j, chain_ix)
                info_i = info_res_first[f"{index_i}{res_i}"]
                info_j = info_res_second[f"{index_j}{res_j}"]
                # Determine whether the residue pair is adjacent
                is_adjacent = self.are_residues_adjacent(index_i, index_j)
                # Adjust atom coordinates based on occupancy
                coord_i = self.adjest_atom_coordinates(is_adjacent, info_i, frame_step)
                coord_j = self.adjest_atom_coordinates(is_adjacent, info_j, frame_step)

                # Pre-filter contacts
                if self.prefilter_contacts(coord_i, coord_j):
                    # Calculate distance between residue pairs
                    dist = self.get_distances(coord_i, coord_j)
                    radius_min = settings['radius_min']
                    radius_max = settings['radius_max']
                    rrcs_score = self.compute_rrcs_jit(dist, radius_max, radius_min)
                else:
                    rrcs_score = 0
                frame_rrcs.append((f"{chain_id}:{index_i}{res_i}", f"{chain_id}:{index_j}{res_j}", rrcs_score))
        return frame_count, frame_rrcs

    @staticmethod
    def check_time_index(begin_index, end_index, step_index):
        """
        Check if the step index is reasonable based on the start and end indices.

        This function is used to validate the time series index step size to ensure it is reasonable,
        and will adjust the step size automatically if it is too large or too small.

        Parameters:
        begin_index: The start index of the time series.
        end_index: The end index of the time series.
        step_index: The step size of the time series index.

        Returns:
        step_index: The adjusted step size.
        """
        index_multiple = int((end_index - begin_index) / step_index) + 1
        if index_multiple <= 0:
            logging.error("The end time is earlier than the start time")
            sys.exit()
        elif index_multiple <= STEP_INDEX_LOWER_LIMIT:
            step_index = int((end_index - begin_index) / 5)
            logging.warning(
                f'Step size too large, automatically reducing step size to {step_index}.')
        elif index_multiple > STEP_INDEX_WARN_LIMIT:
            logging.warning(
                f'The number of frames to be computed exceeds {STEP_INDEX_WARN_LIMIT}. ' 
                + 'While we will not modify the step size, this will require a long computation time.')
        else:
            logging.info(f'Step size is {step_index}, Computation steps are {index_multiple}.')
        
        if begin_index < 0:
            begin_index = 0
            logging.warning('Initial index is less than 0, reset to 0.')
        if step_index <= 0:
            step_index = 1
            logging.warning('Step size is less than 1, reset to 1.')

        return begin_index, step_index

    @timing_decorator
    def analyze_contacts(self, basic_settings):
        """
        Analyze residue-residue contacts in molecular dynamics simulation trajectories.
        
        Parameters:
        basic_settings: dict
            Configuration dictionary containing simulation trajectory and analysis parameters.
            
        Returns:
        all_frame_rrcs: dict
            Dictionary of residue-residue contact information for all frames.
        """
        # Start timing
        global_start = timeit.default_timer()

        # Initialize dictionary to store residue-residue contact information for all frames
        all_frame_rrcs = {}

        # Obtain the molecular dynamics simulation trajectory
        md_traj = basic_settings['md_traj']
        protein = basic_settings['protein']


        # Load the residues of interest based on the residue pair configuration
        member_first, member_second = self.load_residues(basic_settings['res_pairs'])

        info_first = self.get_residue_info(protein, basic_settings['traj_chains'], member_first)
        info_second = self.get_residue_info(protein, basic_settings['traj_chains'], member_second)

        # Use serial processing when the number of CPUs is less than 1
            # Analyze a single frame
        frame_count, frame_rrcs = self.analyze_frame(
                0, 
                info_first, 
                info_second, 
                basic_settings, 
                md_traj)
        all_frame_rrcs[frame_count] = frame_rrcs
        return all_frame_rrcs


def run_pipline():
    """
    The entry point of the program, responsible for executing the analysis and plotting process.
    
    This function first parses the configuration file to obtain basic configuration information,
    then initializes the universe, analyzes the contact data, and finally performs data visualization.
    """
    # Parse the configuration file and get the basic configuration information
    config_parser = ConfigParser()
    basic_config = config_parser.parse_arguments()



    # Initialize the universe based on the basic configuration, and perform initial checks
    initializer = UniverseInitializer(basic_config)
    initializer.run()

    # Create an analyzer instance to analyze the contact data in the simulation
    analyzer = RRCSAnalyzer()
    all_frame_rrcs = analyzer.analyze_contacts(initializer.basic)

    # print(all_frame_rrcs)

    # Create a data visualizer instance and perform data visualization
    DataVisualizer(initializer.basic, all_frame_rrcs).run()


def main():
    """
    The main function of the program, responsible for executing the entire analysis and visualization process.
    """
    global_start = timeit.default_timer()
    run_pipline()
    elapsed = timeit.default_timer() - global_start
    logging.info(
        "The GMX_RRCS program has finished running, with a total elapsed time of "
        + colored(f"{elapsed:.6f} ", 'green')
        + "seconds."
        )


if __name__ == "__main__":

    main()
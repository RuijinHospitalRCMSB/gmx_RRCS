# -*- coding: utf-8 -*-

"""
The `utilities.py` module integrates essential logging, performance profiling, 
and exception handling features aimed at enhancing code robustness and efficiency.

- Logging: Offers `log_error` and `log_warning` functions for recording error and 
  warning messages respectively, color-coded for easy issue identification.
  
- Performance Profiling: The `timing_decorator` decorator automatically measures 
  and logs function execution times, aiding in the optimization of performance 
  bottlenecks.
  
- Exception Handling: Defines custom exceptions `InputFileError`, 
  `ParameterWrongError`, and `ResidueIndexError` for precise capture and handling 
  of specific error scenarios.

This module streamlines error management, boosts self-diagnostic capabilities of 
programs, and assists in optimizing code execution through performance monitoring.
"""


import logging
import timeit
import numpy as np
from numba import jit

from functools import wraps
from termcolor import colored

def log_error(error_type: str, message: str):
    """Log and print the error message."""
    logging.error(message)
    print(colored(f"Error: {error_type}", "red"), message)


def log_warning(warning_type: str, message: str):
    """Log and print the warning message."""
    logging.warning(message)
    print(colored(f"Warning: {warning_type}", "yellow"), message)


def is_alnum_space(s):
    """
    Check if the string s consists only of alphanumeric characters and spaces.

    :param s: str - The string to check.
    :return: bool - True if the string contains only letters, digits, spaces and '-'; False otherwise.
    """
    return all(c.isalnum() or c.isspace() or (c=='-') for c in s)


def timing_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Define an inner function to call the original function with its arguments
        def wrapped_func():
            return func(*args, **kwargs)
        
        # Measure the execution time of the wrapped function using timeit.timeit
        elapsed_time = timeit.timeit(wrapped_func, number=1)
        logging.info(
            f"Function '{func.__name__}' took "
            + colored(f"{elapsed_time:.6f} ", "green")
            + "seconds to complete."
            )
        return wrapped_func()
    return wrapper


def print_nstep_time(n_steps, initial_time, step=5000):
    """
    Prints the elapsed time at specified calculation steps.

    This function is designed to be called at regular intervals during a lengthy calculation process, to log the current
    calculation progress and the time taken. It helps in understanding the computational efficiency and estimating the
    remaining time.

    Parameters:
    - n_steps: int, the current step number of the calculation process.
    - initial_time: float, the start time of the calculation process, usually obtained by calling timeit.default_timer().
    - step: int, the frequency of logging the time, i.e., the function is called to log the time every step steps.
    """
    # Log the progress every {n_steps} frames
    if n_steps % step == 0:
        # Calculate the elapsed time since the start of the calculation
        elapsed = timeit.default_timer() - initial_time
        logging.info(f"Calculated to frame {n_steps}, took {elapsed:.2f} seconds")


class InputFileError(FileNotFoundError):
    pass


class ResidueIndexError(Exception):
    def __init__(self):
        super().__init__(
        colored(f"Error: ", "red"),
        + "The ResidueIndex file format is incorrect."
    )
        
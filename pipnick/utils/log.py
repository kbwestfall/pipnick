from pathlib import Path
import importlib.resources
import logging
import json
from IPython import embed

import numpy as np

def default_log_name(class_name):
    """
    Construct the default log name provided the script class.

    Parameters
    ----------
    class_name : str
        Typically class.__name__ for the calling class.

    Returns
    -------
    str
        Name of the log file.
    """
    return f"{class_name.split('.')[-1]}.log"


def adjust_global_logger(log_level='INFO', log_file=None):
    """
    Adjusts the global logging configuration to use a specific log level and
    output file.

    Parameters
    ----------
    log_level : int, str, optional

        String or integer representation of logging level to set for the
        console handler (default is 'INFO').  Integer mapping is:
        1='CRITICAL', 2='ERROR', 3='WARNING', 4='INFO', 5='DEBUG'.

    log_file : Path, str, optional
        Name for the log file.  If None, no file is written.
    """
    if isinstance(log_level, str):
        _log_level = log_level 
    else:
        log_levels = {1:'CRITICAL', 2:'ERROR', 3:'WARNING', 4:'INFO', 5:'DEBUG'}
        _log_level = log_levels[np.clip(log_level, 1, 5)]

    log_config_file =  importlib.resources.files('pipnick.utils') \
        / ('logging_config_termonly.json' if log_file is None else 'logging_config.json')

    # Adjust the configuration for the file name and log level
    with log_config_file.open('rb') as f:
        config = json.load(f)
    config['handlers']['console']['level'] = _log_level

    if log_file is not None:
        _log_file = Path(log_file).absolute()
        config['handlers']['file']['filename'] = str(_log_file)

    # Configure logging with the loaded configuration
    logging.config.dictConfig(config)


def log_astropy_table(table):
    """
    Formats an Astropy Table into a string for logging purposes.

    Parameters
    ----------
    table : astropy.table.Table
        The Astropy Table to be formatted.

    Returns
    -------
    str
        The formatted string representation of the table.
    """
    # Convert each line of the table's formatted output into a string
    print_table = ""
    for line in table.pformat_all():
        print_table += f"{line}\n"
    return print_table

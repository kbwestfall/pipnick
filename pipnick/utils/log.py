def log_astropy_table(table):
    """
    Formats an Astropy Table into a string for logging purposes.

    Parameters
    ----------
    table : :class:`astropy.table.Table`
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

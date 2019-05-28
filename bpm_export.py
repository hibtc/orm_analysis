#! /usr/bin/env python
"""
Script to read BPM profile exports in .CSV format, calculate mean and standard
error using our own function, and then print out all entries in a simple TXT
format.

Usage:
    ./bpm_export.py <CSV_FILES>...
"""

from docopt import docopt

import pandas as pd
from time import strptime, mktime


def weighted_mean(x, w):
    """Calculate weighted mean of a vector x with weights w."""
    return sum(x * w) / sum(w)


def csv_to_tab(fnames):
    """Read a list of CSV files, calculate mean and standard error from the
    beam profiles using the given function, and return as pandas DataFrame."""
    dicts = []
    for fname in fnames:
        header, custom, data = parse_csv_export(fname)
        time = mktime(strptime(header['Startzeit'], '%d.%m.%Y %H:%M:%S'))
        dicts.append({
            'bpm': header['Gerät'],
            'x_full': weighted_mean(data['Drahtposition x'], data['X']),
            'y_full': weighted_mean(data['Drahtposition y'], data['Y']),
            'x_fwhm': custom['Schwerpunkt X'],
            'y_fwhm': custom['Schwerpunkt Y'],
            'cycle': header['Zyklus ID'],
            'time': round(time),
        })
    return pd.DataFrame(dicts, columns=[
        'bpm',
        'x_full', 'y_full',
        'x_fwhm', 'y_fwhm',
        'cycle', 'time',
    ])


def parse_csv_export(fname):
    """Parses a single CSV file, returns ``(header, custom, data)`` with
    header and custom given as dicts, and data as pandas.DataFrame."""
    with open(fname, encoding='latin1') as f:
        header = dict(_parse_section(f, '<HEADER>', '</HEADER>'))
        custom = dict(_parse_section(f, '<CUSTOM>', '</CUSTOM>'))
        data = pd.read_csv(f, sep=';')
        data.set_index('Drahtnummer')
        return header, custom, data


def _parse_section(lines, start, end):
    """Internal function, parses a <START></END> delimited section in a .CSV
    file, iterates over the items."""
    for line in lines:
        if line.strip() == start:
            break
    for line in lines:
        if line.strip() == end:
            break
        yield line.strip().split(';')


def main(args=None):
    opts = docopt(__doc__, args)
    df = csv_to_tab(opts['<CSV_FILES>'])
    print(df.to_string())


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

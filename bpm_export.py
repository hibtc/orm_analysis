"""
Script to read BPM profile exports in .CSV format, calculate mean and standard
error using our own function, and then print out all entries in a simple TXT
format.

Usage:
    ./bpm_export.py <CSV_FILES>...
"""

import pandas as pd
from time import strptime, mktime


def weighted_mean(x, w):
    """Calculate weighted mean and standard error according to
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Bootstrapping_validation
    """
    w_tot = sum(w)
    x_mean = sum(x * w) / w_tot
    sig_sq = sum(w**2 * (x - x_mean)**2) / w_tot * len(x) / (len(x) - 1)
    return x_mean, sig_sq**0.5


def csv_to_tab(fnames, mean_func=weighted_mean):
    """Read a list of CSV files, calculate mean and standard error from the
    beam profiles using the given function, and return as pandas DataFrame."""
    dicts = []
    for fname in fnames:
        header, custom, data = parse_csv_export(fname)
        x, sigx = mean_func(data['Drahtposition x'], data['X'])
        y, sigy = mean_func(data['Drahtposition y'], data['Y'])
        time = mktime(strptime(header['Startzeit'], '%d.%m.%Y %H:%M:%S'))
        dicts.append({
            'bpm': header['Gerät'],
            'x': x,
            'y': y,
            'sigx': sigx,
            'sigy': sigy,
            'cycle': header['Zyklus ID'],
            'time': round(time),
        })
    return pd.DataFrame(
        dicts, columns=['bpm', 'x', 'y', 'sigx', 'sigy', 'cycle', 'time'])


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


if __name__ == '__main__':
    import sys
    df = csv_to_tab(sys.argv[1:])
    print(df.to_string())

#! /usr/bin/env python
"""
Script to read BPM profile exports in .CSV format, calculate mean and standard
error using our own function, and then print out all entries in a simple TXT
format.

Usage:
    ./bpm_export.py [-m METHOD] <CSV_FILES>...

Options:
    -m METHOD, --mean METHOD        Set method to calculate mean value and
                                    standard error. Possible values are
                                    {full, fwhm} [default: weighted]
"""

from docopt import docopt

import pandas as pd
from time import strptime, mktime
from math import sqrt, log


# conversion factor from FWHM to standard deviation for normal distributions:
fwhm_to_stddev = 2 * sqrt(2 * log(2))


def builtin_mean_func(data, custom):
    return {
        'x': custom['Schwerpunkt X'],
        'y': custom['Schwerpunkt Y'],
        'sigx': custom['FWHM X'] * fwhm_to_stddev,
        'sigy': custom['FWHM Y'] * fwhm_to_stddev,
        'errx': 0,
        'erry': 0,
    }


def weighted_mean_func(data, custom):
    x, sigx, errx = weighted_mean(data['Drahtposition x'], data['X'])
    y, sigy, erry = weighted_mean(data['Drahtposition y'], data['Y'])
    return {'x': x, 'sigx': sigx, 'errx': errx,
            'y': y, 'sigy': sigy, 'erry': erry}


def weighted_mean(x, w):
    """Calculate weighted mean and standard error according to
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Bootstrapping_validation
    """
    w_tot = sum(w)
    x_mean = sum(x * w) / w_tot
    stddev = sum(w * (x - x_mean)**2)**0.5 / w_tot
    sig_sq = sum(w**2 * (x - x_mean)**2) / w_tot * len(x) / (len(x) - 1)
    return x_mean, stddev, sig_sq**0.5


def csv_to_tab(fnames, mean_func=weighted_mean):
    """Read a list of CSV files, calculate mean and standard error from the
    beam profiles using the given function, and return as pandas DataFrame."""
    dicts = []
    for fname in fnames:
        header, custom, data = parse_csv_export(fname)
        time = mktime(strptime(header['Startzeit'], '%d.%m.%Y %H:%M:%S'))
        dicts.append({
            'bpm': header['Gerät'],
            'cycle': header['Zyklus ID'],
            'time': round(time),
            **mean_func(data, custom)
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


MEAN_FUNCS = {
    'full': weighted_mean_func,
    'fwhm': builtin_mean_func,
}


def main(args=None):
    opts = docopt(__doc__, args)
    df = csv_to_tab(
        opts['<CSV_FILES>'],
        MEAN_FUNCS[opts['--method']])
    print(df.to_string())


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

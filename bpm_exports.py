#! /usr/bin/env python
"""
Script to read BPM profile exports in .CSV format, calculate mean and standard
error using our own function, and then print out all entries in a simple TXT
format.

Usage:
    ./bpm_exports.py export <CSV_FILES>...
    ./bpm_exports.py diffhist [-o FILE] <CSV_FILES>...

Options:
    -o FILE, --output FILE          Set image output file
"""

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

from time import strptime, mktime
from math import floor, ceil


def weighted_mean(x, w, rcond=8):
    """Calculate weighted mean of a vector x with weights w."""
    return sum(x * w) / sum(w)


def normal_distribution(x, amplitude, sigma, mu=0, c=0):
    return amplitude * np.exp(-(x - mu)**2 / (2 * sigma**2)) + c


def normal_mean(x, w):
    w_tot = sum(w)
    x_mean = sum(x * w) / w_tot
    x_var = abs(sum(w * (x - x_mean)**2) / w_tot)
    init_vals = [max(w), x_var**0.5, x_mean, min(w)]
    best_vals, covar = curve_fit(
        normal_distribution, x, w, p0=init_vals)
    amplitude, sigma, mu, c = best_vals
    return mu



def csv_to_tab(fnames):
    """Read a list of CSV files, calculate mean and standard error from the
    beam profiles using the given function, and return as pandas DataFrame."""
    dicts = []
    for fname in fnames:
        header, custom, data = parse_csv_export(fname)
        fwhm_x = float(custom['FWHM X'])
        fwhm_y = float(custom['FWHM Y'])
        if fwhm_x == -9999 or fwhm_y == -9999:
            continue
        time = mktime(strptime(header['Startzeit'], '%d.%m.%Y %H:%M:%S'))
        dicts.append({
            'bpm': header['Gerät'],
            'x_full': weighted_mean(data['Drahtposition x'], data['X']),
            'y_full': weighted_mean(data['Drahtposition y'], data['Y']),
            'x_fwhm': float(custom['Schwerpunkt X']),
            'y_fwhm': float(custom['Schwerpunkt Y']),
        #   'x_gauss': normal_mean(data['Drahtposition x'], data['X']),
        #   'y_gauss': normal_mean(data['Drahtposition y'], data['Y']),
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


def diff_histogram(df, bpm, save_to=None):
    """Show a histogram of ``mean_func(x) - mean_acs(x)``, and estimate
    systematic/statistical error."""
    import matplotlib.pyplot as plt
    fig = plt.figure()
    fig.suptitle(bpm)

    for i, axis_name in enumerate('xy'):
        xdata = df[axis_name + '_fwhm']
        ydata = df[axis_name + '_full']

        xmin = floor(min(xdata.min(), ydata.min()))
        xmax = ceil(max(xdata.max(), ydata.max()))

        ax = fig.add_subplot(1, 2, i + 1)
        ax.set_title(axis_name)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(xmin, xmax)
        ax.set_aspect(1)
        ax.plot(df[axis_name + '_fwhm'], df[axis_name + '_full'], 'o')
        ax.plot([xmin, xmax], [xmin, xmax])

    if save_to is None:
        plt.show(fig)
    else:
        fig.savefig(save_to, dpi=200)


def main(args):
    opts = docopt(__doc__, args)
    df = csv_to_tab(opts['<CSV_FILES>'])

    if opts['export']:
        print(df.to_string())

    elif opts['diffhist']:
        for bpm in df.bpm.unique():
            diff_histogram(
                df[df.bpm == bpm], bpm, save_to=opts['--output'])


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

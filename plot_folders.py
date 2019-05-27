#! /usr/bin/env python
"""
Generate plots from previously computed orbit/response data.

Usage:
    ./plot_folders.py [-b BASE] [-m MEASURED] FOLDERS...

Options:
    -m FOLDER, --measured FOLDER    Folder with measured orbits
    -b FOLDER, --base FOLDER        Folder with base model orbits (no errors)
"""

import os

import matplotlib as mpl
import matplotlib.figure
import numpy as np

from docopt import docopt


def main(args=None):
    mpl.rc('text', usetex=False)
    opts = docopt(__doc__, args)
    measured = read_folder(opts['--measured'])
    base = read_folder(opts['--base'])
    for folder in opts['FOLDERS']:
        process_folder(folder, measured or {}, base or {})


def process_folder(folder, measured, base):
    for fname, data in read_folder(folder).items():
        print(fname)
        title = os.path.split(folder)[1] + '/' + fname
        fig = plot_orbit(title, data, measured.get(fname), base.get(fname))
        fig.tight_layout()
        fig.savefig(os.path.join(folder, fname + '.png'))


def read_folder(folder):
    return {
        normalize_name(name): load_table(os.path.join(folder, name))
        for name in os.listdir(folder)
        if name.endswith('.delta') or name.endswith('.orbit')
    } if folder else {}


def normalize_name(fname):
    basename, ext = os.path.splitext(fname)
    return basename[:-2] + ext if basename.endswith('-0') else fname


def load_table(fname):
    return np.genfromtxt(fname)[:, 1:]


def plot_orbit(title, orbit, measured=None, base=None):
    fig = mpl.figure.Figure()
    fig.suptitle(title)
    for i, d in enumerate("xy"):
        ax = fig.add_subplot(1, 2, i + 1)
        ax.set_xlabel(r"monitor position [m]")
        if d == 'x':
            ax.set_ylabel(r"orbit $x$ [mm]")
        else:
            ax.yaxis.tick_right()
        if measured is not None:
            ax.errorbar(
                measured[:, 0],
                measured[:, 1 + 2*i],
                measured[:, 2 + 2*i],
                label=d + " measured")
        ax.plot(
            orbit[:, 0],
            orbit[:, i + 1],
            label=d + " model")
        if base is not None:
            ax.plot(
                base[:, 0],
                base[:, i + 1],
                label=d + " base model")
        ax.legend()
    return fig


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

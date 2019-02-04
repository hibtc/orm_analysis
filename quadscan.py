
import glob
import itertools
import yaml

import numpy as np
import matplotlib.pyplot as plt


def mean_values(data):
    """
    Return NxMx2 array of averaged monitor measurements for the
    N optics, M monitors and 2 spatial dimensions x/y.
    """
    monitors = data['monitors']
    return np.array([
        np.mean([
            [
                [readout[monitor]['posx'],
                 readout[monitor]['posy']]
                for monitor in monitors
            ]
            for record in group
            for readout in [record['readout']]
        ][2:], axis=0)
        for optics, group in itertools.groupby(
                data['records'], key=lambda r: r['optics'])
    ])


def stddevs(data):
    """
    Return NxMx2 array of standard deviations of monitor measurements
    for the N optics, M monitors and 2 spatial dimensions x/y.
    """
    monitors = data['monitors']
    return np.array([
        np.var([
            [
                [readout[monitor]['posx'],
                 readout[monitor]['posy']]
                for monitor in monitors
            ]
            for record in group
            for readout in [record['readout']]
        ][0:], axis=0, ddof=0) ** 0.5
        for optics, group in itertools.groupby(
                data['records'], key=lambda r: r['optics'])
    ])


def read_file(filename):
    """Read a text file, return as string."""
    with open(filename) as f:
        return f.read()


def main():
    filenames = glob.glob('data/2019-01-20_quadscan/*/*/*_X.yml')
    raw_data = [
        yaml.safe_load(read_file(filename))
        for filename in filenames
    ]

    means = [mean_values(data) for data in raw_data]
    errors = [stddevs(data) for data in raw_data]

    kl_values = [
        np.array([
            sum(optics.values())
            for optics, group in itertools.groupby(
                    data['records'], key=lambda r: r['optics'])
        ])
        for data in raw_data
    ]

    for kl, mean, err in zip(kl_values, means, errors):
        plt.errorbar(kl, mean[:, 0, 0], err[:, 0, 0])
        plt.show()


if __name__ == '__main__':
    main()

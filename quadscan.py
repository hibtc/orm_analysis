
import glob
import itertools
import yaml

import numpy as np
import matplotlib.pyplot as plt


def mean_values(data):
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
        for value, group in itertools.groupby(
                data['records'], key=lambda r: r['optics'])
    ])


def stddevs(data):
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
        ][1:], axis=0, ddof=1) ** 0.5
        for value, group in itertools.groupby(
                data['records'], key=lambda r: r['optics'])
    ])


def read_file(filename):
    with open(filename) as f:
        return f.read()


def main():
    filenames = glob.glob('2019-01-20_quadscan/*/*/*_X.yml')
    raw_data = [
        yaml.safe_load(read_file(filename))
        for filename in filenames
    ]

    means = [mean_values(data) for data in raw_data]
    errors = [stddevs(data) for data in raw_data]

    values = [
        np.array([
            sum(value.values())
            for value, group in itertools.groupby(
                    data['records'], key=lambda r: r['optics'])
        ])
        for data in raw_data
    ]

    for vals, mean, err in zip(values, means, errors):
        plt.errorbar(vals, mean[:, 0, 0], err[:, 0, 0])
        plt.show()


if __name__ == '__main__':
    main()

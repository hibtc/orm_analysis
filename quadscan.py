#! /usr/bin/env python

import glob
import itertools
import yaml

import numpy as np
import matplotlib.pyplot as plt

from cpymad.madx import Madx


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
    filenames = glob.glob(
        '../data/orm/2019-01-20_quadscan/M8-E108-F1-I9-G1/*/*_X.yml')
    raw_data = [
        yaml.safe_load(read_file(filename))
        for filename in filenames
    ]

    m = Madx(stdout=False)
    m.verbose()
    m.call('../hit_models/hht3/sequence.madx', chdir=True)
    m.command.beam()
    m.use('hht3')
    #m.call('../hit_models/hht3/strengths0.madx', chdir=True)

    for data in raw_data:

        mean = mean_values(data)
        err = stddevs(data)
        kl = np.array([
            sum(optics.values())
            for optics, group in itertools.groupby(
                    data['records'], key=lambda r: r['optics'])
        ])

        mon = data['monitors'][0]
        knob, = data['optics'][0].keys()
        quad = [elem.name for elem in m.sequence.hht3.expanded_elements
                if elem.base_name == 'quadrupole'
                and knob in m.expr_vars(elem.defs.k1)][0]

        m.globals.update(data['base_optics'])

        for i, x in enumerate('xy'):
            #plt.subplot(2, 1, i+1)
            plt.title(f"pos{x}_{mon}({knob})")
            plt.ylabel("{x} [m]")
            plt.xlabel(f"{knob} [$m^{{-1}}$]")

            plt.errorbar(kl, mean[:, 0, i], err[:, 0, i], label=x)

            for pos in np.linspace(-0.002, 0.002, 5):
                modelled = np.array([
                    track(m, optics, range=f'{quad}/{mon}', **{x: pos})[i]
                    for optics, group in itertools.groupby(
                            data['records'], key=lambda r: r['optics'])
                ])
                i0 = np.argmin(kl)
                modelled += mean[i0, 0, i] - modelled[i0]

                plt.plot(kl, modelled, label=f'${x}_0$={pos}')

            plt.legend()
            plt.show()


def track(madx, optics, **kwargs):
    madx.globals.update(optics)
    tw = madx.twiss(betx=1, bety=1, **kwargs)
    return [tw.x[-1], tw.y[-1]]


if __name__ == '__main__':
    main()

"""
Convert ORM measurement files to numpy/gnuplot compatible files.

Usage:
    ./measured_orbits.py [-m MODEL] [-o DIR] [RECORDS]...

Options:
    -m MODEL, --model MODEL     path to model file
    -o DIR, --output DIR        output base directory

Arguments:
    RECORDS                     path to recorded measurement files
"""

import os

from docopt import docopt
import numpy as np

from madgui.model.orm import Analysis
from util import format_table, format_strengths


def main(args=None):
    opts = docopt(__doc__, args)

    model_file = opts['--model'] or '../hit_models/hht3'
    record_files = (
        opts['RECORDS'] or
        ['2018-10-20-orm_measurements/M8-E108-F1-I9-G1'])

    if len(record_files) == 1 and os.path.isdir(record_files[0]):
        default_prefix = record_files[0] + '_/'
        record_files = record_files[0] + '/*.yml'
    else:
        default_prefix = 'orbits'
    prefix = (opts['--output'] or default_prefix) + '/'

    os.makedirs(prefix, exist_ok=True)

    with Analysis.app(model_file, record_files) as ana:
        elements = ana.model.elements
        records = ana.measured.records
        strengths = ana.measured.strengths
        orbits = {}
        for (monitor, knob), (strength, orbit, variance) in records.items():
            orbits.setdefault(knob, {}) \
                  .setdefault(strength, {})[monitor] = (orbit, variance ** 0.5)

        nan = np.ones((2, 2)) * np.nan
        base_orbit = np.array([
            list(orbits[None][None].get(monitor, nan))
            for monitor in ana.monitors
        ])

        for i, knob in enumerate([None] + ana.knobs):
            by_strength = orbits[knob]
            deltas = sorted(by_strength)
            for j, strength in enumerate(deltas):
                by_monitor = by_strength[strength]
                names = ('name', 's/m', 'x/mm', 'x_err/mm', 'y/mm', 'y_err/mm')
                align = 'lrrrrr'
                formats = [''] + 5 * ['9.5f']
                orbit_table = np.array([
                    list(by_monitor.get(monitor, nan))
                    for monitor in ana.monitors
                ])

                text = format_table(names, align, formats, [
                    (monitor, elements[monitor].position, *values.T.flat)
                    for monitor, values in zip(ana.monitors, orbit_table * 1e3)
                ])
                basename = (f'{prefix}{i}_base' if knob is None else
                            f'{prefix}{i}_{knob}-{j}')
                with open(f'{basename}.orbit', 'wt') as f:
                    f.write(text)

                str_data = format_strengths(
                    strengths if knob is None else {
                        knob: strength,
                        knob + '__orig': strengths[knob],
                        knob + '__delta': strength - strengths[knob],
                    })
                with open(f'{basename}.str', 'wt') as f:
                    f.write(str_data)

                if knob is not None:
                    response = np.dstack((
                        (orbit_table[:, 0, :] - base_orbit[:, 0, :]),
                        (orbit_table[:, 1, :] ** 2 + base_orbit[:, 1, :] ** 2)
                        ** 0.5,
                    ))

                    text = format_table(names, align, formats, [
                        (monitor, elements[monitor].position, *values.flat)
                        for monitor, values in zip(
                                ana.monitors, response * 1e3)
                    ])
                    with open(f'{basename}.delta', 'wt') as f:
                        f.write(text)


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

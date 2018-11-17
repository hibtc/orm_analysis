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
from itertools import starmap

from docopt import docopt

from madgui.model.orm import Analysis


def main(args=None):
    opts = docopt(__doc__, args)

    model_file = opts['--model'] or '../hit_models/hht3'
    record_files = (
        opts['RECORDS'] or
        ['../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1'])

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

        null = (0, 0), (0, 0)
        for i, knob in enumerate([None] + ana.knobs):
            by_strength = orbits[knob]
            deltas = sorted(orbits[knob])
            for j, strength in enumerate(deltas):
                by_monitor = by_strength[strength]
                names = ('name', 's/m', 'x/mm', 'dx/mm', 'y/mm', 'dy/mm')
                align = 'lrrrrr'
                formats = [''] + 5 * ['.5f']
                orbit_table = format_table(names, align, formats, [
                    (monitor, elements[monitor].position,
                     orbit[0] * 1e3, error[0] * 1e3,
                     orbit[1] * 1e3, error[1] * 1e3)
                    for monitor in ana.monitors
                    for orbit, error in [by_monitor.get(monitor, null)]
                ])
                basename = (f'{prefix}orbit-{i}_base' if knob is None else
                            f'{prefix}orbit-{i}_{knob}-{j}')
                with open(f'{basename}.txt', 'wt') as f:
                    f.write(orbit_table)

                str_data = format_strengths(
                    strengths if knob is None else {
                        knob: strength,
                        knob + '__orig': strengths[knob],
                        knob + '__delta': strength - strengths[knob],
                    })
                with open(f'{basename}.str', 'wt') as f:
                    f.write(str_data)


def format_table(names, align, formats, table, sep='   '):
    header = ['# ' + names[0]] + list(names[1:])
    cells = [header] + [
        list(starmap(format, zip(items, formats)))
        for items in table
    ]
    widths = list(map(max, zip(*[
        map(len, items) for items in cells
    ])))
    return '\n'.join([
        sep.join(starmap(adjust, zip(items, align, widths)))
        for items in cells
    ])


_adjust = {
    'l': str.ljust,
    'r': str.rjust,
}


def adjust(s, a, w):
    return _adjust[a](s, w)


def format_strengths(data):
    return ''.join([
        '{} = {!r};\n'.format(k, v)
        for k, v in data.items()
    ])


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

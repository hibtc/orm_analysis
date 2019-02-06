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

from orm_util import Analysis
from util import format_table, format_strengths


def main(args=None):
    opts = docopt(__doc__, args)

    model_file = opts['--model'] or '../hit_models/hht3'
    record_files = (
        opts['RECORDS'] or
        ['data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1'])

    if len(record_files) == 1 and os.path.isdir(record_files[0]):
        default_prefix = record_files[0] + '_/'
        record_files = record_files[0] + '/*.yml'
    else:
        default_prefix = 'orbits'
    prefix = (opts['--output'] or default_prefix) + '/'

    os.makedirs(prefix, exist_ok=True)

    with Analysis.app(model_file, record_files) as ana:
        elements = ana.model.elements
        strengths = ana.measured.strengths
        orbits = ana.measured.orbits
        errors = ana.measured.stddev
        base_orbit = orbits[:, :, 0]
        base_error = errors[:, :, 0]

        for i, optic in enumerate(ana.optics):
            orbit = np.dstack([
                orbits[:, :, i],
                errors[:, :, i],
            ])
            names = ('name', 's/m', 'x/mm', 'x_err/mm', 'y/mm', 'y_err/mm')
            align = 'lrrrrr'
            formats = [''] + 5 * ['9.5f']

            text = format_table(names, align, formats, [
                (monitor, elements[monitor].position, *values.flat)
                for monitor, values in zip(ana.monitors, orbit * 1e3)
            ])
            label = next(iter(optic))[0] if optic else 'base'
            basename = f'{prefix}{i}_{label}'
            with open(f'{basename}.orbit', 'wt') as f:
                f.write(text)

            str_data = format_strengths(dict(optic or strengths))
            with open(f'{basename}.str', 'wt') as f:
                f.write(str_data)

            if optic:
                response = np.dstack((
                    (orbit[:, :, 0] - base_orbit),
                    (orbit[:, :, 1] ** 2 + base_error ** 2)
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

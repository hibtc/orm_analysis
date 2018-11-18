"""
Compute model orbits and responses.

Usage:
    ./plot_responses.py [-m MODEL] [-o DIR] [-e ERRORS] [RECORDS]...

Options:
    -m MODEL, --model MODEL         path to model file
    -o DIR, --output DIR            output base directory
    -e FILE, --error FILE           file with error definitions

Arguments:
    RECORDS                         path to recorded measurement files
                                    we use the same deltas as given in the
                                    measurements
"""

import os

import numpy as np
from docopt import docopt

from madgui.model.errors import parse_error, apply_errors, Param
from madgui.model.orm import Analysis
import madgui.util.yaml as yaml

from util import format_table


def get_orbit(model, errors, values):
    """Get x, y vectors, with specified errors."""
    madx = model.madx
    madx.command.select(flag='interpolate', clear=True)
    with apply_errors(model, errors, values):
        tw_args = model._get_twiss_args(table='orm_tmp')
        twiss = madx.twiss(**tw_args)
    return np.stack((twiss.x, twiss.y)).T


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

    if opts['--error']:
        spec = yaml.load_file(opts['--error'])
    else:
        spec = {}
    errors = list(map(parse_error, spec.keys()))
    values = list(spec.values())

    with Analysis.app(model_file, record_files) as ana:
        model = ana.model
        strengths = ana.measured.strengths
        monitors = ana.monitors
        knobs = ana.knobs
        elements = model.elements
        deltas = {
            knob: strength - strengths[knob]
            for (_, knob), (strength, __, ___) in ana.measured.records.items()
            if knob is not None
        }

        model.madx.eoption(add=True)
        model.update_globals(strengths.items())

        sel = [model.elements.index(m) for m in monitors]
        orbits = np.array([get_orbit(model, errors, values)] + [
            get_orbit(model, [Param(knob)] + errors, [deltas[knob]] + values)
            for knob in knobs
        ])[:, sel, :]

        for i, (knob, orbit) in enumerate(zip([None] + knobs, orbits)):

            names = ('name', 's/m', 'x/mm', 'y/mm')
            align = 'lrrr'
            formats = [''] + 3 * ['9.5f']

            text = format_table(names, align, formats, [
                (monitor, elements[monitor].position, *values.flat)
                for monitor, values in zip(monitors, orbit * 1e3)
            ])
            basename = (f'{prefix}{i}_base' if knob is None else
                        f'{prefix}{i}_{knob}')
            with open(f'{basename}_model.orbit', 'wt') as f:
                f.write(text)

            if knob is not None:
                response = orbit - orbits[0]

                text = format_table(names, align, formats, [
                    (monitor, elements[monitor].position, *values.flat)
                    for monitor, values in zip(monitors, response * 1e3)
                ])
                with open(f'{basename}_model.delta', 'wt') as f:
                    f.write(text)


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

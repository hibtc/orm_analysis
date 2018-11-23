"""
Compute model orbits and responses.

Usage:
    ./model_orbits.py [-m MODEL] [-o DIR] [-i ERRORS] [-e ERRORS] [RECORDS]...

Options:
    -m MODEL, --model MODEL         path to model file
    -o DIR, --output DIR            output base directory
    -i FILE, --init FILE            file with global error definitions
    -e FILE, --error FILE           file with error definitions

Arguments:
    RECORDS                         path to recorded measurement files
                                    we use the same deltas as given in the
                                    measurements
"""

import os
from shutil import copyfile

import numpy as np
from docopt import docopt

from madgui.model.errors import parse_error, apply_errors, Param
from madgui.model.orm import Analysis
import madgui.util.yaml as yaml

from util import format_table, format_strengths


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
        ['2018-10-20-orm_measurements/M8-E108-F1-I9-G1'])

    if len(record_files) == 1 and os.path.isdir(record_files[0]):
        default_prefix = record_files[0] + '_'
        record_files = record_files[0] + '/*.yml'
    else:
        default_prefix = 'orbits'

    if opts['--init']:
        g_spec = yaml.load_file(opts['--init'])
    else:
        g_spec = {}
    g_errors = list(map(parse_error, g_spec.keys()))
    g_values = list(g_spec.values())

    if opts['--error']:
        specs = yaml.load_file(opts['--error'])
        default_prefix += '/' + os.path.splitext(
            os.path.basename(opts['--error']))[0]
    else:
        specs = [{}]
    if isinstance(specs, dict):
        specs = [specs]

    prefix = (opts['--output'] or default_prefix) + '/'
    os.makedirs(prefix, exist_ok=True)
    if opts['--init']:
        copyfile(opts['--init'], f'{prefix}spec_init.yml')
    if opts['--error']:
        copyfile(opts['--error'], f'{prefix}spec_error.yml')

    with Analysis.app(model_file, record_files) as ana:
        model = ana.model
        strengths = ana.measured.strengths

        model.madx.eoption(add=True)
        model.update_globals(strengths.items())

        for i, spec in enumerate(specs):
            errors = list(map(parse_error, spec.keys())) + g_errors
            values = list(spec.values()) + g_values
            errname = '_' + repr(errors[0]) if len(spec) == 1 else ''
            output_orbits(ana, f'{prefix}/model_{i}{errname}/', errors, values)


def output_orbits(ana, prefix, errors, values):
    os.makedirs(prefix, exist_ok=True)

    spec = format_strengths(dict(zip(map(repr, errors), values)))
    with open(f'{prefix}errors.txt', 'wt') as f:
        f.write(spec.replace('=', ':'))

    model = ana.model
    strengths = ana.measured.strengths
    monitors = ana.monitors
    knobs = ana.knobs
    elements = model.elements
    records = ana.measured.records
    deltas = {
        knob: strength - strengths[knob]
        for (monitor, knob), (strength, orbit, error) in records.items()
        if knob is not None
    }

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
        with open(f'{basename}.orbit', 'wt') as f:
            f.write(text)

        if knob is not None:
            response = orbit - orbits[0]

            text = format_table(names, align, formats, [
                (monitor, elements[monitor].position, *values.flat)
                for monitor, values in zip(monitors, response * 1e3)
            ])
            with open(f'{basename}.delta', 'wt') as f:
                f.write(text)


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

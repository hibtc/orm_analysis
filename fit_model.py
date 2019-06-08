#! /usr/bin/env python
"""
Fit model to measurements using the given error parameters.

Usage:
    ./fit_model.py [-m MODEL] (-e ERRORS) [-b] [-o FOLDER] [<MEASURED>...]

Options:
    -b, --backtrack                 Use backtracking
    -m MODEL, --model MODEL         Model path
    -e ERRORS, --errors ERRORS      File containing initial values for errors
    -o FOLDER, --output FOLDER      Output folder [default: results]

Arguments:
    <MEASURED>                      YAML files with madgui measurements
"""

import os
import sys
from scipy.optimize import Bounds
from docopt import docopt

import madgui.util.yaml as yaml

from orm_util import Analysis, parse_errors


def main(args=None):

    opts = docopt(__doc__, args)
    prefix = opts['--output']
    os.makedirs(prefix, exist_ok=True)

    record_files = (
        opts['<MEASURED>'] or
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

    model_path = opts['--model'] or '../hit_models/hht3'
    ana = Analysis.app(model_path, record_files)

    if opts['--backtrack']:
        from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
        ana.ensure_monitors_available(from_monitors)
        ana.setup_backtracking(ana.extrapolate(from_monitors, to='#e'))

    ana.init()

    ana.plot_orbit(save_to=prefix+'/0-init')

    errors = parse_errors(yaml.load_file(opts['--errors']))

    options = dict(
        algorithm='lstsq',
        mode='xy',
        iterations=8,
        delta=1e-3,
        bounds=Bounds(-0.001, 0.001),
        fourier=False,
        rcond=1e-6,
    )

    result = ana.fit(errors, save_to=prefix+'/2-fit.txt', **options)

    print(result)

    ana.plot_monitors(save_to=prefix+'/3-final')
    ana.plot_orbit(save_to=prefix+'/3-final')
    ana.plot_steerers(save_to=prefix+'/3-final')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

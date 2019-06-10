#! /usr/bin/env python
"""
Fit model to measurements using the given error parameters.

Usage:
    ./fit_model.py [-m MODEL] [-b] [-o FOLDER] [--stderr VAL]
                   [-absolute | --hybrid] [-i ERRORS]
                   <ERRORS> [<MEASURED>...]

Options:
    -b, --backtrack                 Use backtracking
    -m MODEL, --model MODEL         Model path
    -o FOLDER, --output FOLDER      Output folder [default: results]
    --absolute                      Fit absolute orbits
    --hybrid                        Fit absolute orbit + orbit response
    --stderr VALUE                  Set BPM position stderr [default: 0.0003]
    -i ERRORS                       Set initial/global errors

Arguments:
    <ERRORS>                        YAML file containing initial error values
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

    BPM_ERR = float(opts['--stderr'])

    record_files = (
        opts['<MEASURED>'] or
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

    model_path = opts['--model'] or '../hit_models/hht3'
    ana = Analysis.app(model_path, record_files)
    ana.model.update_twiss_args(dict(x=0, y=0, px=0, py=0))
    if opts['--absolute']:
        ana.mode = 'abs'
    elif opts['--hybrid']:
        ana.mode = 'hybrid'
    else:
        ana.mode = 'orm'

    ana.measured.stderr = (ana.measured.stderr**2 + BPM_ERR**2)**0.5

    if opts['-i']:
        init_errors = parse_errors(yaml.load_file(opts['-i']))
        ana.apply_errors(init_errors.keys(), init_errors.values())

    if opts['--backtrack']:
        from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
        ana.ensure_monitors_available(from_monitors)
        ana.setup_backtracking(ana.extrapolate(from_monitors, to='#e'))

    ana.init()

    errors = parse_errors(yaml.load_file(opts['<ERRORS>']))

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


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

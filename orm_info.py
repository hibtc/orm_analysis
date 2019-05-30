#! /usr/bin/env python
"""
orm_show

Usage:
    ./orm_info.py [-o DIR] [-m MODEL] [-b] [--base-orbit]
                  [-e ERRORS]... [RECORDS]...

Options:
    -e FILE, --errors FILE      apply errors from file
    -o DIR, --output DIR        set output directory [default: results]
    -m PATH, --model PATH       model path [default: ../hit_models/hht3]
    -b, --backtrack             use backtracking
    --base-orbit                plot unmodified orbit as well

Arguments:
    MODEL                       madx model file
    RECORDS                     measured ORM records
"""

import os

from docopt import docopt

import madgui.util.yaml as yaml
from madgui.model.errors import import_errors

from orm_util import Analysis


def main(args=None):
    opts = docopt(__doc__, args)

    prefix = opts['--output'] + '/'
    os.makedirs(prefix, exist_ok=True)

    model = opts['--model']
    record_files = (
        opts['RECORDS'] or
        'data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

    ana = Analysis.app(model, record_files)

    if opts['--base-orbit']:
        ana.init()
        base_orm = ana.model_orbits
    else:
        base_orm = None

    for filename in opts['--errors']:
        import_errors(ana.model, yaml.load_file(filename))

    ana.init({})

    if opts['--backtrack']:
        ana._init_twiss = ana.extrapolate(
            ['t3dg2g', 't3dg1g', 't3df1'], to='#s')

    ana.info()
    ana.plot_monitors(base_orm=base_orm, save_to=f'{prefix}mon_response')
    ana.plot_steerers(base_orm=base_orm, save_to=f'{prefix}mon_response')
    ana.plot_orbit(save_to=f'{prefix}orbit')


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))

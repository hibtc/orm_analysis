#! /usr/bin/env python
"""
orm_show

Usage:
    ./orm_info.py [-e ERRORS]... [-o DIR] [MODEL] [RECORDS]...

Options:
    -e FILE, --errors FILE      apply errors from file
    -o DIR, --output DIR        set output directory [default: results]

Arguments:
    MODEL                       madx model file
    RECORDS                     measured ORM records
"""

import os

from docopt import docopt

import madgui.util.yaml as yaml
from madgui.model.errors import import_errors

from orm_util import Analysis

opts = docopt(__doc__)

model = opts['MODEL'] or 'hit_models/hht3'
record_files = (
    opts['RECORDS'] or
    'data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')


with Analysis.app(model, record_files) as ana:

    ana.init()
    #ana.backtrack(['h1dg1g', 'h1dg2g', 'h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g'])
    base_orm = ana.model_orm

    for filename in opts['--errors']:
        import_errors(ana.model, yaml.load_file(filename))
    ana.init({})
    ana.backtrack(['h1dg1g', 'h1dg2g', 'h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g'])
    #ana.backtrack(['h3dg3g', 'b3dg2g', 'b3dg3g'])

    #print(ana.backtrack(['h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g']))

    prefix = opts['--output'] + '/'
    os.makedirs(prefix, exist_ok=True)

    ana.info()
    ana.plot_monitors(['g3dg5g'], base_orm=base_orm)
    #ana.plot_steerers(['h1ms4v'], base_orm=base_orm)
    ana.plot_orbit()
    #ana.plot_monitors(['g3dg3g'], save_to=prefix, base_orm=base_orm)

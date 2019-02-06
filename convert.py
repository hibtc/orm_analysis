"""
Used to convert a YAML file created from the madgui dialog at:

    Online control -> Beam diagnostic -> Offsets -> Calibrate

to the format as required by ``fit_model.py``, i.e. the format that is created
by the dialog at:

    Online control -> ORM measurement
"""

import sys
import os

import yaml

from orm_util import _convert_orm_export


def convert(filename):
    with open(filename) as f:
        data = yaml.safe_load(f)
    os.makedirs('conv/' + os.path.split(filename)[0], exist_ok=True)
    converted = _convert_orm_export(data)
    text = yaml.safe_dump(converted)
    with open(f'conv/{filename}', 'wt') as f:
        f.write(text)


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        convert(fname)

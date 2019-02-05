"""
Used to convert a YAML file created from the madgui dialog at:

    Online control -> Beam diagnostic -> Offsets -> Calibrate

to the format as required by ``fit_model.py``, i.e. the format that is created
by the dialog at:

    Online control -> ORM measurement
"""

import sys
import os
import itertools

import yaml


def convert(filename):
    with open(filename) as f:
        data = yaml.safe_load(f)

    monitors = data['monitors']
    knobs = {
        knob
        for optic in data['optics']
        for knob in optic
    }

    converted_data = yaml.safe_dump({
        'model': data['base_optics'],
        'monitors': monitors,
        'knobs': list(knobs),
        'sequence': 'hht3',
        'records': [
            {
                'optics': optic,
                'shots': [
                    {monitor: [r['posx'], r['posy'], r['envx'], r['envy']]
                     for monitor in monitors
                     for r in [shot['readout'][monitor]]}
                    for shot in group
                ],
            }
            for optic, group in itertools.groupby(
                data['records'], key=lambda r: r['optics'])
        ]
    })


    os.makedirs('conv/' + os.path.split(filename)[0], exist_ok=True)
    with open(f'conv/{filename}', 'wt') as f:
        f.write(converted_data)


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        convert(fname)

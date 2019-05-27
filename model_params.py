"""
Generate a list of error parameters for the given model.

Usage:
    ./model_params.py [-m MODEL] [-a ATTR]... [-b BLACKLIST]... 
                      [--twiss] [--ealign] [--knobs]

Options:
    --ealign                        List alignment errors
    --knobs                         List knob errors
    --twiss                         List twiss args errors
    -m MODEL, --model MODEL         Path to model
    -a ATTR, --attr ATTR            List element attribute (e.g. "sbend->e1/e2")
    -b ELEM, --blacklist ELEM       Blacklist element

Example:
    ./model_params.py -m hht3 \
        -a quadrupole->k1 \
        -a sbend->angle/e1/e2/k1 \
        -a translation->x/y/px/py \
        -b s4mu1e -b s4mu2e \
        -b s3me2e -b s3mu1a -b s3ms1v \
        -b s0qg4d -b s0qg1f

Missing:
    - efcomp (field errors)
"""

from madgui.core.app import init_app
from madgui.core.session import Session
from madgui.core.config import load as load_config

from docopt import docopt


def main(args=None):
    opts = docopt(__doc__, args)
    model_path = opts['--model'] or '../hit_models/hht3'

    init_app(['madgui'])

    blacklist = opts['--blacklist'] or [
        's4mu1e', 's4mu2e',
        's3me2e', 's3mu1a', 's3ms1v',
        's0qg4d', 's0qg1f',
    ]

    def blacklisted(name):
        return name in blacklist or any(b in name for b in blacklist)

    config = load_config(isolated=True)
    session = Session(config)
    session.load_model(model_path, stdout=False)
    model = session.model()

    by_type = get_elements_by_type(model)

    if opts['--twiss']:
        print('')
        print('############')
        print('## TWISS  ##')
        print('############')
        print('\n'.join(get_twiss_args_errors()))

    if opts['--ealign']:
        print('')
        print('############')
        print('## EALIGN ##')
        print('############')
        for base_name in sorted(by_type):
            print('')
            print('# {}'.format(base_name))
            for elem in by_type[base_name]:
                if not blacklisted(elem):
                    print('\n'.join(get_ealign_errors(elem)))

    if opts['--knobs']:
        print('')
        print('############')
        print('## KNOBS  ##')
        print('############')

        for base_name in sorted(by_type):
            print('')
            print('# {}'.format(base_name))
            for elem in by_type[base_name]:
                if not blacklisted(elem):
                    print('\n'.join(get_knob_errors(model, elem)))

    if opts['--attr']:
        print('')
        print('############')
        print('## ATTRS  ##')
        print('############')
        for spec in opts['--attr']:
            type_, attrs = spec.split('->', 1)
            attrs = attrs.split('/')
            for elem in model.elements:
                if elem.base_name == type_:
                    print('\n'.join(
                        '{}->{}'.format(elem.name, attr) for attr in attrs))


def get_elements_by_type(model):
    supported_types = {
        'sbend', 'quadrupole', 'hkicker', 'vkicker', 'kicker',
        'solenoid', 'multipole', 'srotation',
    }
    by_type = {}
    for elem in model.elements:
        if elem.base_name in supported_types:
            by_type.setdefault(elem.base_name, []).append(elem)
    return by_type


def get_twiss_args_errors():
    return ['x', 'y', 'px', 'py']


def get_knob_errors(model, elem):
    yield from map('δ{}'.format, model.get_elem_knobs(elem))


def get_ealign_errors(elem):
    return ['{}<{}>'.format(elem.name, attr)
            for attr in ('dx', 'dy', 'ds', 'dphi', 'dtheta', 'dpsi')]


if __name__ == '__main__':
    import sys; sys.exit(main(sys.argv[1:]))

"""
Generate a list of error parameters for the given model.

Currently prints:

    - initial particle coordinates
    - ealign errors
    - knobs (prints scaling errors)

Missing:
    - efcomp (field errors)
"""

from madgui.core.app import init_app
from madgui.core.session import Session
from madgui.core.config import load as load_config


def get_grouped_params(model):
    supported_elements = {
        'sbend', 'quadrupole', 'hkicker', 'vkicker', 'kicker',
        'solenoid', 'multipole', 'srotation'
    }
    by_type = {}
    for elem in model.elements:
        if elem.base_name in supported_elements:
            by_type.setdefault(elem.base_name, []).append(elem)

    yield ''
    yield '############'
    yield '## TWISS  ##'
    yield '############'
    yield 'x'
    yield 'y'
    yield 'px'
    yield 'py'

    yield ''
    yield '############'
    yield '## EALIGN ##'
    yield '############'
    for base_name in sorted(by_type):
        yield ''
        yield '# {}'.format(base_name)
        for elem in by_type[base_name]:
            yield from ealign(elem)

    yield ''
    yield '############'
    yield '## KNOBS  ##'
    yield '############'
    for base_name in sorted(by_type):
        yield ''
        yield '# {}'.format(base_name)
        for elem in by_type[base_name]:
            yield from map('δ{}'.format, model.get_elem_knobs(elem))


def ealign(elem):
    for attr in ('dx', 'dy', 'ds', 'dphi', 'dtheta', 'dpsi'):
        yield '{}<{}>'.format(elem.name, attr)


def main(model_file='hit_models/hht3'):
    init_app(['madgui'])

    config = load_config(isolated=True)
    with Session(config) as session:
        session.load_model(
            model_file,
            stdout=False)
        model = session.model()
        print('\n'.join(get_grouped_params(model)))


if __name__ == '__main__':
    import sys; sys.exit(main(*sys.argv[1:]))

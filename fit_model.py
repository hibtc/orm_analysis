import sys
from madgui.model.errors import parse_error
from scipy.optimize import Bounds

from cpymad.types import VAR_TYPE_DIRECT

from orm_util import Analysis


def parse_errors(names):
    return list(map(parse_error, names))


record_files = (
    sys.argv[1:] or
    '../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

with Analysis.app('../hit_models/hht3', record_files) as ana:

    ana.model = ana.model.reversed()

    ana.init()
    ana.backtrack(['h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g'])

    madx = ana.model.madx
    elems = ana.model.elements
    bends = [elem for elem in elems if elem.base_name == 'sbend']
    quads = [elem for elem in elems if elem.base_name == 'quadrupole']
    kicks = [elem for elem in elems if elem.base_name.endswith('kicker')]
    patch = [elem for elem in elems if elem.base_name == 'translation']

    e_orbit = parse_errors(['x', 'y', 'px', 'py'])

    e_ealign = parse_errors([
        f'{elem.name}<{err}>'
        for elem in quads + bends
        for err in ['dx', 'dy', 'ds', 'dphi', 'dtheta', 'dpsi']
    ])

    e_quad_k1 = parse_errors([f'δ{elem.name}->k1' for elem in quads])
    e_bend_k1 = parse_errors([f'Δ{elem.name}->k1' for elem in bends])
    e_bend_angle = parse_errors([f'Δ{elem.name}->angle' for elem in bends])

    e_kick = parse_errors([
        f'δ{knob}'
        for knob, par in ana.model.globals.cmdpar.items()
        if knob.split('_')[0] in ('ax', 'ay', 'dax')
        and par.var_type == VAR_TYPE_DIRECT
    ])

    e_patch = parse_errors([
        f'{elem.name}->{err}'
        for elem in patch
        for err in ('x', 'y', 'px', 'py')
    ])

    errors = sum((
        e_orbit,
        e_ealign,
        e_quad_k1,
        e_bend_k1,
        e_bend_angle,
        e_kick,
        e_patch,
    ), [])
    monitors = ana.monitors

    options = dict(
        algorithm='svd',
        mode='xy',
        iterations=30,
        delta=1e-5,
        bounds=Bounds(-0.001, 0.001),
        fourier=True,
    )

    ana.fit(errors, monitors, save_to='plots/2-fit.txt', **options)

    ana.plot_monitors(monitors, save_to='plots/3-final')
    ana.plot_orbit(save_to='plots/3-final')

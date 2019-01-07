import sys
from madgui.model.errors import parse_error, apply_errors, Param
from madgui.online.orbit import fit_particle_readouts, Readout
from scipy.optimize import Bounds
import numpy as np

from cpymad.types import VAR_TYPE_DIRECT

from orm_util import Analysis


def parse_errors(names):
    return list(map(parse_error, names))


record_files = (
    sys.argv[1:] or
    '../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')


def extrapolate_orbit(measured, i_knob, model, from_monitors, to='#e'):
    # TODO: in the more general case, we would also need to set the strengths
    # corresponding to i_knob
    return fit_particle_readouts(model, [
        Readout(monitor, *measured.orm[index, :, i_knob])
        for monitor in from_monitors
        for index in [measured.monitors.index(monitor.lower())]
    ], to=to)[0][0]


with Analysis.app('../hit_models/hht3', record_files) as ana:

    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
    final_orbits = {
        knob: extrapolate_orbit(ana.measured, i, ana.model, from_monitors)
        for i, knob in enumerate([None] + ana.knobs)
    }
    reverse_init_orbits = {
        knob: {'x': -orbit['x'], 'px': orbit['px'],
               'y': orbit['y'], 'py': -orbit['py']}
        for knob, orbit in final_orbits.items()
    }

    def get_orbit(errs, vals, knob):
        model = ana.model
        madx = model.madx
        madx.command.select(flag='interpolate', clear=True)

        deltas = ana.deltas
        errors = [Param(knob)] + errs if knob else errs
        values = [deltas[knob]] + vals if knob else vals

        with apply_errors(model, errors, values):
            tf = madx.twiss(
                table='orm_tmp', sequence='hht3',
                betx=1, bety=1, **reverse_init_orbits[knob])
            return np.stack((tf.x, tf.y)).T

    ana._get_orbit = get_orbit

    ana.model.reverse()
    ana.model.update_twiss_args(reverse_init_orbits[None])
    ana.measured.orm[:, 0, :] *= -1

    ana.init()

    madx = ana.model.madx
    elems = ana.model.elements
    bends = [elem for elem in elems if elem.base_name == 'sbend']
    quads = [elem for elem in elems if elem.base_name == 'quadrupole']
    kicks = [elem for elem in elems if elem.base_name.endswith('kicker')]
    patch = [elem for elem in elems if elem.base_name == 'translation']

    e_ealign = parse_errors([
        f'{elem.name}<{err}>'
        for elem in quads + bends
        for err in ['dx', 'dy', 'ds', 'dphi', 'dtheta', 'dpsi']
    ])

    e_quad_k1 = parse_errors([f'δ{elem.name}->k1' for elem in quads])
    e_bend_k1 = parse_errors([f'Δ{elem.name}->k1' for elem in bends])
    e_bend_angle = parse_errors([f'Δ{elem.name}->angle' for elem in bends])

    e_bend_angle = parse_errors([
        f'Δ{knob}'
        for knob, par in ana.model.globals.cmdpar.items()
        if knob.split('_')[0] == 'dax'
        and par.var_type == VAR_TYPE_DIRECT
    ])

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

    e_patch = parse_errors([
        'x_g3tx1',
        'y_g3tx1',
        'px_g3tx1',
        'py_g3tx1',
    ])

    # TODO:
    # - fit iteratively on last element using differences
    # - add final offset in the end
    # - impose further constraint = same end position/angle

    # TODO: plot all orbits

    errors = sum((
    #   e_ealign,
    #   e_quad_k1,
    #   e_bend_k1,
        e_bend_angle,
    #   e_kick,
        e_patch,
    ), [])

    errors = parse_errors([
        'Δx_g3tx1', 'Δy_g3tx1',
        'Δpx_g3tx1', 'Δpy_g3tx1',
        'Δdax_g3mu1', 'Δdax_g3mu2', 'Δdax_g3mu3',
        'Δdax_b3mu1', 'Δdax_b3mu2',
    ]) + e_kick

    monitors = ana.monitors

    options = dict(
        algorithm='lstsq',
        mode='xy',
        iterations=30,
        delta=1e-5,
        bounds=Bounds(-0.001, 0.001),
        fourier=True,
    )

    ana.fit(errors, monitors, save_to='plots/2-fit.txt', **options)

    ana.plot_monitors(monitors, save_to='plots/3-final')
    ana.plot_orbit(save_to='plots/3-final')

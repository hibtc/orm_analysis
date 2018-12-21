import sys
from madgui.model.errors import parse_error, apply_reverse_errors, Param
from madgui.model.orm import Analysis
from madgui.online.orbit import fit_particle_readouts, Readout
from scipy.optimize import Bounds
import numpy as np


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
    ])


with Analysis.app('../hit_models/hht3', record_files) as ana:

    ana.init()

    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
    final_orbits = {
        knob: extrapolate_orbit(ana.measured, i, ana.model, from_monitors)
        for i, knob in enumerate([None] + ana.knobs)
    }

    # prepare the reversed sequence
    ana.model.backtrack(x=0, y=0)

    def get_orbit(errs, vals, knob):
        model = ana.model
        madx = model.madx
        madx.command.select(flag='interpolate', clear=True)

        deltas = ana.deltas
        errors = [Param(knob)] + errs if knob else errs
        values = [deltas[knob]] + vals if knob else vals

        # TODO reverse errors
        with apply_reverse_errors(model, errors, values):
            ti = final_orbits[knob]
            tf = madx.twiss(
                table='orm_tmp', sequence=model.backseq,
                betx=1, bety=1,
                x=-ti['x'], px=ti['px'],
                y=ti['y'], py=-ti['py'])
            return np.stack((-tf.x, tf.y)).T[::-1]

    ana._get_orbit = get_orbit

    madx = ana.model.madx
    elems = ana.model.elements

    orbit_errors = parse_errors([
        'chr->x', 'chr->y', 'chr->px', 'chr->py',
    ])

    ealign_errors = parse_errors([
        '{}<{}>'.format(elem.name, err)
        for elem in elems
        if elem.base_name in ('quadrupole', 'sbend')
        for err in ('dx', 'dy', 'ds')
        #for err in ('dx', 'dy', 'ds', 'dphi', 'dtheta', 'dpsi')
    ])

    quad_errors = parse_errors([
        'δ{}->k1'.format(elem.name)
        for elem in elems
        if elem.base_name in ('quadrupole', 'sbend')
    ])

    axgeo_errors = parse_errors([
        'Δ{}'.format(elem.defs.angle)
        for elem in elems
        if elem.base_name == 'sbend'
        and isinstance(elem.defs.angle, str)
        and elem.defs.angle.startswith('axgeo_')
    ])

    bend_errors = parse_errors([
        'δ{}'.format(knob)
        for elem in elems
        if elem.base_name == 'sbend'
        for knob in (set(madx.expr_vars(elem.defs.k0)) -
                     set(madx.expr_vars(elem.defs.angle)))
    ])

    kick_errors = parse_errors([
        'δ{}->kick'.format(elem.name)
        for elem in elems
        if elem.base_name in ('hkicker', 'vkicker')
    ])

    errors = (
        #orbit_errors +
        ealign_errors +
        quad_errors +
        axgeo_errors +
        bend_errors +
        kick_errors
    )

    # Nelder-Mead
    # CG
    # BFGS
    # L-BFGS-B
    # TNC
    options = dict(
        algorithm='svd',
        #method='Nelder-Mead',
        mode='xy',
        iterations=30,
        delta=1e-5,
        bounds=Bounds(-0.001, 0.001),
        fourier=True,
    )

    ana.fit(errors, ana.monitors, save_to='plots/2-fit.txt', **options)

    ana.plot_monitors(['g3dg5g'], save_to='plots/3-final')
    ana.plot_orbit(save_to='plots/3-final')

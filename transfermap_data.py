import os
import sys
from madgui.online.orbit import fit_particle_readouts, Readout
import numpy as np

import matplotlib.pyplot as plt

from orm_util import Analysis

record_files = (
    sys.argv[1:] or
    'data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')


def extrapolate_orbit(measured, i_knob, model, from_monitors, to='#e'):
    # TODO: in the more general case, we would also need to set the strengths
    # corresponding to i_knob
    return fit_particle_readouts(model, [
        Readout(monitor, *measured.orbits[index, :, i_knob])
        for monitor in from_monitors
        for index in [measured.monitors.index(monitor.lower())]
    ], to=to)[0][0]


with Analysis.app('../hit_models/hht3', record_files) as ana:

    exclude_knobs = ['ay_g3mw1', 'ax_g3mw2', 'dax_g3mu3']

    obs_el = 'g3dg5g'
    obs_idx = ana.monitors.index(obs_el)

    # estimate particle coordinates at ISO center:
    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
    final_orbits = {
        knob: extrapolate_orbit(ana.measured, i, ana.model, from_monitors)
        for i, knob in enumerate([None] + ana.knobs)
    }
    reverse_init_orbits = {
        knob: (-orbit['x'], orbit['px'],
               orbit['y'], -orbit['py'])
        for knob, orbit in final_orbits.items()
    }

    ana.model.update_globals(ana.measured.strengths)
    ana.model.reverse()
    ana.model.update_twiss_args(dict(zip(
        ('x', 'px', 'y', 'py'), reverse_init_orbits[None])))
    T_model = ana.model.sectormap('#s', obs_el)
    S = T_model[[0, 2]]

    x_iso = np.array([
        reverse_init_orbits[knob]
        for i, knob in enumerate([None] + ana.knobs)
        #if knob is None or knob.lower() not in exclude_knobs
    ])

    y_obs = np.array([
        ana.measured.orbits[obs_idx, :, i]
        for i, knob in enumerate([None] + ana.knobs)
        #if knob is None or knob.lower() not in exclude_knobs
    ])
    y_obs[:, 0] *= -1       # -x

    def twiss(x, px, y, py):
        tw = ana.model.track_one(x=x, px=px, y=y, py=py, range=f'#s/{obs_el}')
        return [tw.x[-1], tw.y[-1]]

    y_twiss = np.array([
        twiss(*reverse_init_orbits[knob])
        for knob in [None] + ana.knobs
        #if knob is None or knob.lower() not in exclude_knobs
    ]).T

    np.savetxt(
        'results/transfer_particles.txt', np.hstack((x_iso, y_obs)),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g")

    # fit transfermap
    N = len(x_iso)
    H = np.ones((N, 1))
    x_iso = np.hstack((x_iso, H))
    result_x = np.linalg.lstsq(x_iso, y_obs[:, 0], rcond=1e-8)
    result_y = np.linalg.lstsq(x_iso, y_obs[:, 1], rcond=1e-8)
    tx = result_x[0]
    ty = result_y[0]
    T = np.vstack((tx, ty))

    TT = np.hstack((S[:, 0:4], S[:, [6]]))

    np.savetxt('results/transfer_map.txt', np.vstack((T, TT)))

    y_fit = T @ x_iso.T

    patch = (np.linalg.lstsq(H, (y_obs.T - y_twiss)[0], rcond=1e-8)[0],
             np.linalg.lstsq(H, (y_obs.T - y_twiss)[1], rcond=1e-8)[0])
    print(patch)

    # model

    #madx = ana.model.madx
    #madx.input('''
    #    select, flag=makethin, class=sbend, slice=50, thick=false;
    #    makethin, sequence=hht3, style=teapot, makedipedge=true;
    #    use, sequence=hht3;
    #''')

    #def track(x, px, y, py):
    #    madx.input(f'''
    #        track, onepass=true, onetable=true;
    #        start, x={x}, y={y}, px={px}, py={py};
    #        observe, place={obs_el};
    #        run, turns=1;
    #        endtrack;
    #    ''')
    #    tab = madx.table.trackone
    #    return [tab.x[1], tab.y[1]]

    #y_track = np.array([
    #    track(*reverse_init_orbits[knob])
    #    for knob in [None] + ana.knobs
    #    if knob is None or knob.lower() not in exclude_knobs
    #]).T

    y_model = y_twiss

    resid_e = y_fit - y_obs.T
    resid_m = y_model - y_obs.T

    def red_chisq(x, err=1, ddof=0):
        res = (x / err).flatten()
        return np.dot(res, res) / (len(res) - ddof)

    print(red_chisq(resid_e, 0.5e-3, 8))
    print(red_chisq(resid_m, 0.5e-3, 8))

    # plot

    os.makedirs('results/tm_plots', exist_ok=True)

    for iy, y_ax in enumerate('xy'):
        fig = plt.figure()
        for ix, x_ax in enumerate(('x', 'px', 'y', 'py')):
            ax = fig.add_subplot(2, 2, 1 + ix)
            ax.set_title(f'{y_ax}({x_ax})')
            I = np.argsort(x_iso[:, ix])
            ax.plot(x_iso[I, ix], y_obs[I, iy], label=f'measured')
            ax.plot(x_iso[I, ix], y_model[iy, I], label=f'model')
            ax.plot(x_iso[I, ix], y_fit[iy, I], label='fit')
            ax.plot(x_iso[I, ix], y_model[iy, I] + patch[iy], label=f'patch')
            #ax.plot(x_iso[I, ix], y_twiss[iy, I], 'o-', label='twiss')
            #ax.plot(x_iso[I, ix], y_track[iy, I], '-', label='track')
        legend = ax.legend(loc='upper center', fancybox=True,
                           bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
        fig.savefig(f'results/tm_plots/{y_ax}.png', bbox_inches='tight')
        plt.clf()

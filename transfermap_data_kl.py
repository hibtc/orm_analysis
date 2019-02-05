import os
import sys
from madgui.online.orbit import fit_particle_readouts, Readout
import numpy as np

import matplotlib.pyplot as plt

from orm_util import Analysis


record_files = (
    sys.argv[1:] or
    './data/2019-01-20_quadscan/M8-E108-F1-I9-G1/0-*/*.yml')


def extrapolate_orbit(measured, i_optic, model, from_monitors, to='#e'):
    # TODO: in the more general case, we would also need to set the strengths
    # corresponding to i_optic
    return fit_particle_readouts(model, [
        Readout(monitor, *measured.orbits[index, :, i_optic])
        for monitor in from_monitors
        for index in [measured.monitors.index(monitor.lower())]
    ], to=to)[0][0]


with Analysis.app('../hit_models/hht3', record_files) as ana:

    exclude_knobs = ['ay_g3mw1', 'ax_g3mw2', 'dax_g3mu3']

    obs_el = 'g3dg5g'
    obs_idx = ana.monitors.index(obs_el)

    # estimate particle coordinates at ISO center:
    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']
    final_orbits = [
        extrapolate_orbit(ana.measured, i, ana.model, from_monitors)
        for i, optic in enumerate(ana.optics)
    ]
    reverse_init_orbits = [
        (-orbit['x'], orbit['px'], orbit['y'], -orbit['py'])
        for orbit in final_orbits
    ]

    ana.model.update_globals(ana.measured.strengths)
    ana.model.reverse()
    ana.model.update_twiss_args(dict(zip(
        ('x', 'px', 'y', 'py'), reverse_init_orbits[0])))
    T_model = ana.model.sectormap('#s', obs_el)
    S = T_model[[0, 2]]

    x_iso = np.array(reverse_init_orbits)
    y_obs = np.array([
        ana.measured.orbits[obs_idx, :, i]
        for i, _ in enumerate(ana.optics)
        #if knob is None or knob.lower() not in exclude_knobs
    ])
    y_obs[:, 0] *= -1       # -x

    def twiss(x, px, y, py):
        tw = ana.model.track_one(x=x, px=px, y=y, py=py, range=f'#s/{obs_el}')
        return [tw.x[-1], tw.y[-1]]

    y_twiss = np.array([
        twiss(*orbit)
        for orbit in reverse_init_orbits
        #if knob is None or knob.lower() not in exclude_knobs
    ]).T

    base_kl = ana.strengths['kl_g3qd42']
    kl = np.column_vector([
        optic.get('kl_g3qd42', base_kl)
        for optic in ana.optics
    ])

    np.savetxt(
        'results/transfer_particles.txt', np.hstack((x_iso, y_obs, kl)),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g kl_g3qd42")

    # fit transfermap
#   N = len(x_iso)
#   H = np.ones((N, 1))
#   x_iso = np.hstack((x_iso, H))
#   result_x = np.linalg.lstsq(x_iso, y_obs[:, 0], rcond=1e-8)
#   result_y = np.linalg.lstsq(x_iso, y_obs[:, 1], rcond=1e-8)
#   tx = result_x[0]
#   ty = result_y[0]
#   T = np.vstack((tx, ty))

#   TT = np.hstack((S[:, 0:4], S[:, [6]]))

#   np.savetxt('results/transfer_map.txt', np.vstack((T, TT)))

#   y_fit = T @ x_iso.T

#   patch = (np.linalg.lstsq(H, (y_obs.T - y_twiss)[0], rcond=1e-8)[0],
#            np.linalg.lstsq(H, (y_obs.T - y_twiss)[1], rcond=1e-8)[0])
#   print(patch)

    # model

    y_model = y_twiss

#   resid_e = y_fit - y_obs.T
#   resid_m = y_model - y_obs.T

#   def red_chisq(x, err=1, ddof=0):
#       res = (x / err).flatten()
#       return np.dot(res, res) / (len(res) - ddof)

#   print(red_chisq(resid_e, 0.5e-3, 8))
#   print(red_chisq(resid_m, 0.5e-3, 8))

    # plot

    os.makedirs('results/tm_plots_kl', exist_ok=True)

    # TODO: plot _kL value on x axis

    for iy, y_ax in enumerate('xy'):
        fig = plt.figure()
        for ix, x_ax in enumerate(('x', 'px', 'y', 'py')):
            ax = fig.add_subplot(2, 2, 1 + ix)
            ax.set_title(f'{y_ax}({x_ax})')
            I = np.argsort(kl)
            ax.plot(kl[I], y_obs[I, iy], label=f'measured')
            ax.plot(kl[I], y_model[iy, I], label=f'model')
       #    ax.plot(kl[I], y_fit[iy, I], label='fit')
       #    ax.plot(kl[I], y_model[iy, I] + patch[iy], label=f'patch')
            #ax.plot(x_iso[I, ix], y_twiss[iy, I], 'o-', label='twiss')
            #ax.plot(x_iso[I, ix], y_track[iy, I], '-', label='track')
        legend = ax.legend(loc='upper center', fancybox=True,
                           bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
        fig.savefig(f'results/tm_plots/{y_ax}.png', bbox_inches='tight')
        plt.clf()

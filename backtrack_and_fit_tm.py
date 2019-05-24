"""
I wrote this script to visualize and investigate the discrepancies between
model and measurements of G3MU3 (the final 90° gantry dipole).

Usage:
    ./transfermap_data.py [<MEASURED>]

Example:
    ./transfermap_data.py 2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml

It calculates beam momentum at T3DF1 from position measurements at the three
MWPCs and then uses these as initial coordinates to backtrack to G3DG5G. The
reverse sectormap T[T3DF1->G3DG5G] is calculated once using the MAD-X model
and once fitted using a least squares fit.

Finally, it plots the measured available BPM readout (one for each optic) at
G3DG5G versus the beam coordinates at T3DF1, compared to the backtracked
positions.

Output will be stored under `results/tm`:

    x.png                           X at G3DG5G versus X/PX/Y/PY at T3DF1
    y.png                           Y at G3DG5G versus X/PX/Y/PY at T3DF1
    tm_model.txt                    2x5 transfermap (x,px,y,py,1) → (x,y)
    tm_fit.txt                      fitted transfermap
    measured_beam_coordinates.txt   Initial coordinates for backtracking for
                                    each of the available optics
"""

import os
import sys
from madgui.online.orbit import fit_particle_readouts, Readout
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from orm_util import Analysis


def main(record_files):
    record_files = (
        record_files or
        'data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

    ana = Analysis.app('../hit_models/hht3', record_files)

    prefix = 'results/tm/'
    os.makedirs(prefix, exist_ok=True)

    # show results at this element:
    obs_el = 'g3dg5g'
    obs_idx = ana.monitors.index(obs_el)

    # use these bpms to guess the initial conditions for the backtracking run:
    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']

    # pick optics for which we have measured all involved BPMs:
    involved_bpms = [obs_idx] + [ana.monitors.index(m) for m in from_monitors]
    usable_optics = [
        i_optic for i_optic in range(len(ana.optics))
        if not np.any(np.isnan(
            ana.measured.orbits[involved_bpms][:, :, i_optic]))
    ]
    num_optics = len(usable_optics)

    # estimate particle coordinates at ISO center for each of the optics:
    final_orbits = [
        extrapolate_orbit(ana.measured, i_optic, ana.model, from_monitors)
        for i_optic in usable_optics
    ]

    # prepare initial coordinates for backtracking:
    reverse_init_orbits = pd.DataFrame([
        (-orbit['x'], orbit['px'],
         orbit['y'], -orbit['py'])
        for orbit in final_orbits
    ], columns=['x', 'px', 'y', 'py'])

    # measured positions with computed momenta at T3DF1:
    x_iso = np.array(reverse_init_orbits)
    # measured positions at G3DG5G:
    y_obs = ana.measured.orbits[obs_idx][:, usable_optics].T
    y_obs[:, 0] *= -1       # X[reverse] = -X[forward]

    # save measured beam coordinates:
    np.savetxt(
        prefix+'measured_beam_coordinates.txt', np.hstack((x_iso, y_obs)),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g")

    # compute positions at G3DGG by backtracking from G3DG5G:
    ana.model.update_globals(ana.measured.strengths)
    ana.model.reverse()
    ana.model.update_twiss_args(reverse_init_orbits.iloc[0].to_dict())

    T_model = ana.model.sectormap('#s', obs_el)[[0, 2]][:, [0, 1, 2, 3, 6]]
    y_model = np.array([
        track(ana.model, *orbit, range=f'#s/{obs_el}')
        for orbit in x_iso
    ]).T

    # fit transfermap + kick (x/px/y/py/1 -> x/y)
    H = np.ones((num_optics, 1))
    x_iso = np.hstack((x_iso, H))
    T_fit = np.vstack((
        np.linalg.lstsq(x_iso, y_obs[:, 0], rcond=1e-8)[0],
        np.linalg.lstsq(x_iso, y_obs[:, 1], rcond=1e-8)[0],
    ))

    np.savetxt(prefix+'/tm_model.txt', T_model, header='x px y py k')
    np.savetxt(prefix+'/tm_fit.txt', T_fit, header='x px y py k')

    y_fit = T_fit @ x_iso.T

    offset = np.vstack([
        np.linalg.lstsq(H, (y_obs.T - y_model)[0], rcond=1e-8)[0],
        np.linalg.lstsq(H, (y_obs.T - y_model)[1], rcond=1e-8)[0],
    ])
    y_offset = y_model + offset
    print("fitted x/y offsets:", *offset)

    # model
    resid_e = y_fit - y_obs.T
    resid_m = y_model - y_obs.T
    resid_o = y_offset - y_obs.T

    print("red χ² [init]   =", red_chisq(resid_m, err=0.5e-3, ddof=0))
    print("red χ² [fit]    =", red_chisq(resid_e, err=0.5e-3, ddof=8))
    print("red χ² [offset] =", red_chisq(resid_o, err=0.5e-3, ddof=1))

    # plot
    for iy, y_ax in enumerate('xy'):
        fig = plt.figure(figsize=(8, 8))
        fig.suptitle(
            f"backtracked {y_ax} at {obs_el} versus beam coordinates at ISO center")
        for ix, x_ax in enumerate(('x', 'px', 'y', 'py')):
            ax = fig.add_subplot(2, 2, 1 + ix)
            ax.set_title(f'{y_ax}({x_ax})')
            I = np.argsort(x_iso[:, ix])
            ax.plot(x_iso[I, ix], y_obs[I, iy], label=f'measured')
            # ax.plot(x_iso[I, ix], y_model[iy, I], label=f'model')
            ax.plot(x_iso[I, ix], y_fit[iy, I], label='fit')
            ax.plot(x_iso[I, ix], y_offset[iy, I], label=f'model + {offset[iy][0]}')
        ax.legend(loc='upper center', fancybox=True,
                  bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
        fig.savefig(prefix+f'{y_ax}.png', bbox_inches='tight')
        plt.clf()


def extrapolate_orbit(measured, i_optic, model, from_monitors, to='#e'):
    """Extrapolate particle position/momentum from the position measurements
    of the given BPMs ``from_monitors``.

    This function does NOT update optics and is therefore only elligible for
    pure DRIFT sections."""
    # TODO: in the more general case, we would also need to set the strengths
    # corresponding to i_optic
    return fit_particle_readouts(model, [
        Readout(monitor, *measured.orbits[index, :, i_optic])
        for monitor in from_monitors
        for index in [measured.monitors.index(monitor.lower())]
    ], to=to)[0][0]


def track(model, x, px, y, py, range='#s/#e'):
    """Return final (x, y) for a particle with the given initial conditions."""
    tw = model.track_one(x=x, px=px, y=y, py=py, range=range)
    return [tw.x[-1], tw.y[-1]]


def red_chisq(x, err=1, ddof=0):
    res = (x / err).flatten()
    return np.dot(res, res) / (len(res) - ddof)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

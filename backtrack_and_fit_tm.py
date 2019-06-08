#! /usr/bin/env python
"""
I wrote this script to visualize and investigate the discrepancies between
model and measurements of G3MU3 (the final 90° gantry dipole).

Usage:
    ./transfermap_data.py [<MEASURED>...]

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
from math import ceil, floor
import numpy as np

import matplotlib.pyplot as plt

from orm_util import Analysis


def main(record_files):
    record_files = (
        record_files or
        #'../data/orm/*/M8-E1*/*.yml')
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

    ana = Analysis.app('../hit_models/hht3', record_files)

    prefix = 'results/tm/'
    os.makedirs(prefix, exist_ok=True)

    # show results at this element:
    obs_el = 'g3dg5g'
    obs_idx = ana.monitors.index(obs_el)

    # use these bpms to guess the initial conditions for the backtracking run:
    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']

    ana.ensure_monitors_available([obs_el] + from_monitors)
    ana.setup_backtracking(ana.extrapolate(from_monitors, to='#e'))

    num_optics = len(ana.optics)

    # measured positions with computed momenta at T3DF1:
    x_iso = np.array([
        [twiss['x'], twiss['px'], twiss['y'], twiss['py']]
        for optic in ana.optics
        for twiss in [ana._init_twiss[optic]]
    ])
    # measured positions at G3DG5G:
    y_obs = ana.measured.orbits[obs_idx]
    y_err = ana.measured.stddev[obs_idx]

    # save measured beam coordinates:
    np.savetxt(
        prefix+'measured_beam_coordinates.txt', np.hstack((x_iso, y_obs.T)),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g")

    # compute positions at G3DGG by backtracking from G3DG5G:
    ana.model.update_globals(ana.measured.strengths)

    T_model = ana.model.sectormap('#s', obs_el)[[0, 2]][:, [0, 1, 2, 3, 6]]
    y_model = ana.compute_model_orbits()[obs_idx]

    # fit transfermap + kick (x/px/y/py/1 -> x/y)
    H = np.ones((num_optics, 1))
    x_iso = np.hstack((x_iso, H))
    T_fit = np.vstack((
        np.linalg.lstsq(x_iso, y_obs[0], rcond=1e-8)[0],
        np.linalg.lstsq(x_iso, y_obs[1], rcond=1e-8)[0],
    ))
    y_fit = T_fit @ x_iso.T

    np.savetxt(prefix+'/tm_model.txt', T_model, header='x px y py k')
    np.savetxt(prefix+'/tm_fit.txt', T_fit, header='x px y py k')

    # fit constant difference:
    offset = np.vstack([
        np.linalg.lstsq(H, (y_obs - y_model)[0], rcond=1e-8)[0],
        np.linalg.lstsq(H, (y_obs - y_model)[1], rcond=1e-8)[0],
    ])
    y_offset = y_model + offset
    print("fitted x/y offsets:", *offset)

    # estimate fit qualities:
    resid_e = y_fit - y_obs
    resid_m = y_model - y_obs
    resid_o = y_offset - y_obs

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
            ax.plot(x_iso[I, ix], y_obs[iy, I], label=f'measured')
            # ax.plot(x_iso[I, ix], y_model[iy, I], label=f'model')
            ax.plot(x_iso[I, ix], y_fit[iy, I], label='fit')
            ax.plot(x_iso[I, ix], y_offset[iy, I], label=f'model + {offset[iy][0]}')
        ax.legend(loc='upper center', fancybox=True,
                  bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
        fig.savefig(prefix+f'{y_ax}.png', bbox_inches='tight')
        plt.clf()

        # with #optic on the x axis:
        fig = plt.figure(figsize=(8, 8))
        off = offset[iy][0] * 1000
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(f"backtracked {y_ax} at {obs_el} versus optic index")
        ax.set_xlabel(f'optic')
        ax.set_ylabel(f'{y_ax} / mm')
        ax.plot(y_obs[iy]*1000, 'o', label=f'measured')
        # ax.plot(y_model[iy], label=f'model')
        ax.plot(y_fit[iy]*1000, label='fit')
        ax.plot(y_offset[iy]*1000, label=fr'$\mathrm{{model}}{off:+.2f}$')

        ax.legend(fancybox=True, shadow=True)
        fig.savefig(prefix+f'{y_ax}_by_optic.png', bbox_inches='tight')
        plt.clf()

    print("plot measured versus expected X/Y")
    fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 1]})
    fig.suptitle('beam position at G3DG5G', y=0.87)

    for iy, ax, name in zip(range(2), axes, 'XY'):
        xdata = y_obs[iy] * 1000
        ydata = y_model[iy] * 1000
        yerr = y_err[iy] * 1000

        offs = -offset[iy, 0]*1000
        yoff = xdata + offs

        I = np.argsort(xdata)
        ax.errorbar(xdata[I], ydata[I], yerr[I],
                    fmt='none', label="model", color=f'C0')

        ax.plot(xdata[I], yoff[I], '-', label=f'measured ${offs:+.1f}$',
                color='C1')

        ax.set_title(name)
        ax.set_xlabel(f'measured {name} [mm]')

        ax.legend()
        ax.set_aspect(1)

    axes[0].set_ylabel(f'model [mm]')
    fig.savefig(prefix+f'beam_pos_g3dg5g.png', bbox_inches='tight', dpi=400)
    plt.clf()


    print("Plotting measured versus model orbit (same plot)")
    fig = plt.figure()
    fig.suptitle('beam position at G3DG5G', y=0.87)
    ax_x = fig.add_subplot(1, 1, 1)
    ax_R = ax_x.twinx()
    ax_y = ax_R.twiny()
    ax_x.set_xlabel(f'measured X [mm]', color='C0')
    ax_y.set_xlabel(f'measured Y [mm]', color='C1')
    ax_x.set_ylabel(f'model X [mm]', color='C0')
    ax_R.set_ylabel(f'model Y [mm]', color='C1')
    ax_x.tick_params(axis='x', labelcolor='C0')
    ax_x.tick_params(axis='y', labelcolor='C0')
    ax_y.tick_params(axis='x', labelcolor='C1')
    ax_R.tick_params(axis='y', labelcolor='C1')
    ax_x.spines['left'].set_color('C0')
    ax_x.spines['bottom'].set_color('C0')
    ax_x.spines['top'].set_color('C1')
    ax_x.spines['right'].set_color('C1')

    xdata = y_obs * 1000
    ydata = y_model * 1000
    yerr = y_err * 1000

    I = np.argsort(xdata[0])
    ax_x.errorbar(xdata[0, I], ydata[0, I], yerr[0, I],
                  fmt='', label='X', color='C0')

    I = np.argsort(xdata[1])
    ax_y.errorbar(xdata[1, I], ydata[1, I], yerr[1, I],
                  fmt='', label='Y', color='C1')

    #ax_x.set_aspect(1)
    #ax_y.set_aspect(1)

    set_lims(ax_x, xdata[0, I], ydata[0, I])
    set_lims(ax_y, xdata[1, I], ydata[1, I])

    fig.savefig(prefix+f'measured_vs_model.png', bbox_inches='tight', dpi=400)
    plt.clf()


def set_lims(ax, xdata, ydata):
    xmin, xmax = floor(min(xdata) - 0.5), ceil(max(xdata) + 0.5)
    ymin, ymax = floor(min(ydata) - 0.5), ceil(max(ydata) + 0.5)

    size = max(xmax - xmin, ymax - ymin)
    xadd = size - (xmax - xmin)
    yadd = size - (ymax - ymin)

    xmin -= xadd / 2
    xmax += xadd / 2
    ymin -= yadd / 2
    ymax += yadd / 2

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def red_chisq(x, err=1, ddof=0):
    res = (x / err).flatten()
    return np.dot(res, res) / (len(res) - ddof)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

#! /usr/bin/env python
import os
import sys
import numpy as np

import matplotlib.pyplot as plt

from orm_util import Analysis


def main(record_files):
    record_files = (
        record_files or
        'data/2019-01-20_quadscan/M8-E108-F1-I9-G1/3-*/*.yml')

    ana = Analysis.app('../hit_models/hht3', record_files)

    prefix = 'results/kl_dependency/'
    os.makedirs(prefix, exist_ok=True)

    # show results at this element:
    obs_el = 'g3dg5g'
    obs_idx = ana.monitors.index(obs_el)

    # use these bpms to guess the initial conditions for the backtracking run:
    from_monitors = ['t3dg2g', 't3dg1g', 't3df1']

    ana.ensure_monitors_available([obs_el] + from_monitors)
    ana.setup_backtracking(ana.extrapolate(from_monitors, to='#e'))

    # measured positions with computed momenta at T3DF1:
    x_iso = np.array([
        [twiss['x'], twiss['px'], twiss['y'], twiss['py']]
        for optic in ana.optics
        for twiss in [ana._init_twiss[optic]]
    ])
    # measured positions at G3DG5G:
    y_obs = ana.measured.orbits[obs_idx].T

    base_kl = ana.measured.strengths['kl_efg_g3qd42']
    kl = np.array([
        dict(optic).get('kl_efg_g3qd42', base_kl)
        for optic in ana.optics
    ])
    print(kl)
    print(ana.optics)

    np.savetxt(
        prefix+'measured_beam_coordinates.txt', np.hstack((x_iso, y_obs, kl[:, None])),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g kl_g3qd42")

    # compute positions at G3DGG by backtracking from G3DG5G:
    ana.model.update_globals(ana.measured.strengths)
    y_model = ana.compute_model_orbits()[obs_idx]

    for iy, y_ax in enumerate('xy'):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(f'{y_ax}(kL)')
        I = np.argsort(kl)
        ax.plot(kl[I], y_obs[I, iy], label=f'measured')
        ax.plot(kl[I], y_model[iy, I], label=f'model')
        ax.legend(loc='upper center', fancybox=True,
                  bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
        fig.savefig(prefix+f'{y_ax}.png', bbox_inches='tight')
        plt.clf()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

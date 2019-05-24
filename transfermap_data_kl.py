import os
import sys
from madgui.online.orbit import fit_particle_readouts, Readout
import pandas as pd
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

    # pick optics for which we have measured all involved BPMs:
    involved_bpms = [obs_idx] + [ana.monitors.index(m) for m in from_monitors]
    usable_optics = [
        i_optic for i_optic in range(len(ana.optics))
        if not np.any(np.isnan(
            ana.measured.orbits[involved_bpms][:, :, i_optic]))
    ]

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

    base_kl = ana.measured.strengths['kl_efg_g3qd42']
    kl = np.array([
        dict(ana.optics[i_optic]).get('kl_efg_g3qd42', base_kl)
        for i_optic in usable_optics
    ])
    print(kl)
    print(ana.optics)

    np.savetxt(
        prefix+'measured_beam_coordinates.txt', np.hstack((x_iso, y_obs, kl[:, None])),
        header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g kl_g3qd42")

    # compute positions at G3DGG by backtracking from G3DG5G:
    ana.model.update_globals(ana.measured.strengths)
    ana.model.reverse()
    ana.model.update_twiss_args(reverse_init_orbits.iloc[0].to_dict())

    y_model = np.array([
        track(ana.model, *orbit, range=f'#s/{obs_el}')
        for orbit in x_iso
    ]).T

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


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

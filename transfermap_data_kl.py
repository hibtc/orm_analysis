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


ana = Analysis.app('../hit_models/hht3', record_files)

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
])
y_obs[:, 0] *= -1       # -x

def twiss(x, px, y, py):
    tw = ana.model.track_one(x=x, px=px, y=y, py=py, range=f'#s/{obs_el}')
    return [tw.x[-1], tw.y[-1]]

y_twiss = np.array([
    twiss(*orbit)
    for orbit in reverse_init_orbits
]).T

base_kl = ana.strengths['kl_g3qd42']
kl = np.column_vector([
    optic.get('kl_g3qd42', base_kl)
    for optic in ana.optics
])

np.savetxt(
    'results/transfer_particles.txt', np.hstack((x_iso, y_obs, kl)),
    header="x_iso px_iso y_iso py_iso x_g3dg5g y_g3dg5g kl_g3qd42")

y_model = y_twiss

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
    legend = ax.legend(loc='upper center', fancybox=True,
                       bbox_to_anchor=(-0.1, -0.2), shadow=True, ncol=4)
    fig.savefig(f'results/tm_plots/{y_ax}.png', bbox_inches='tight')
    plt.clf()

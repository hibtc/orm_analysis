#! /usr/bin/env python
# Plots:


# For all data:
# backtracked model orbit versus monitors for
#   - entire line
#   - zoom in on gantry


# Backtracked positions versus measurements in G3DG5G


# Measured versus modelled orbit resposne matrix:
#   - before gantry
#   - on gantry


# Fitted orbit response matrix:
#   - before gantry
#   - on gantry


# BPM profile, gaussian, non-gaussian 

import shutil
import os
from glob import glob

import numpy as np

from madgui.online.orbit import fit_particle_readouts, Readout
from madgui.model.madx import Model
from madgui.util.yaml import load_file
from orm_plot import create_twiss_figure


# ====
# DATA
# ====

DATA_PREFIX = '../data/'
HEBT_BPMS = ['h1dg1g', 'h1dg2g', 'h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g']
ISOC_BPMS = ['t3dg2g', 't3dg1g', 't3df1']

# Folders with corresponding '.knobs.str' and '.bpms.yml' files with
# measurements for a single optic:
single_optic_measurements_folders = [
    'correct/2018-06-18-correct/',      # the other dataset is a bit nicer
    'correct/2018-07-03-correct/',
#   'correct/2018-07-15-correct/',      # only has MWPCs, not very useful
]

# Folders with ORM measurements. Each of these folders corresponds to one ORM
# measurement, and contains multiple YAML files, each one corresponding to a
# set of monitors:
orm_measurements_folders = [
    'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/',
    'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G4/',
    'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G4-ISO/',
    'orm/2019-01-20-orm_measurements/M3-E108-F1-I9-G10-ISO/',
    'orm/2019-01-20-orm_measurements/M8-E108-F1-I9-G1-ISO/',
    'orm/2019-01-20-quadscan/M8-E108-F1-I9-G1/',
    'orm/2019-05-11-orm_measurements/M8-E108-F1-I9-G1/',
    'orm/2019-05-11-orm_measurements/M8-E108-F1-I9-G1-LIN/',
]

# CSV files with BPM profile exports:
bpm_profile_exports = [
    'emittance/2017-09-26-emittance/4300/*/*.CSV',
    'emittance/2018-08-05-emittance/4300/*/*.CSV',
    'emittance/2018-08-05-emittance/4400/*/*.CSV',
    'emittance/2019-01-19-emittance/4300/*/*.CSV',
    'emittance/2019-01-19-emittance/4400/*/*.CSV',
    'emittance/2019-04-15-emittance/4300/*/*.CSV',
    'emittance/2019-04-26-emittance/4400/*/*.CSV',
    'orm/2019-05-11-orm_measurements/Gitter-Exporte/*/*.CSV',
]

# what about offset calibrations?
offset_calib_data = [
    'offsets/2018-05-06-offset-calibration/',
]


def plot_single_optic_measurements():
    """Create plots for the 'single optic' measurements (i.e. combinations of
    .knobs.str + .bpms.yml exports)."""
    # Currently, we have only gantry measurements of this type:
    model_name = 'hht3'
    bpm_suffix = '.bpms.yml'
    str_suffix = '.knobs.str'
    bpm_files = [
        fname
        for folder in single_optic_measurements_folders
        for fname in glob(DATA_PREFIX + folder + '*' + bpm_suffix)
    ]

    for bpm_file in bpm_files:
        str_file = strip_suffix(bpm_file, bpm_suffix) + str_suffix
        img_base = strip_suffix(bpm_file, bpm_suffix)

        print("Processing", bpm_file)

        model = Model.load_file(
            f'../hit_models/{model_name}.cpymad.yml', stdout=False)
        model.interpolate = 400
        model.call(str_file)

        # dict with {<monitor>: {'x': <X>, 'y': <Y>, ...}}
        measured = load_file(bpm_file)['monitor']
        measured = sort_dict(measured, key=model.elements.index)

        bpm_names = sorted(measured, key=model.elements.index)
        bpm_data = np.array([
            [model.elements[bpm_name].position,
             measured[bpm_name]['x'],
             measured[bpm_name]['y']]
            for bpm_name in bpm_names
        ])

        model_hebt = interpolate(model, measured, HEBT_BPMS)
        model_isoc = interpolate(model, measured, ISOC_BPMS)

        curves = [
            (bpm_data.T, args('o', label='measured', color='C2')),
            (model_hebt, args('-', label='model', color='C0')),
            (model_isoc, args('-', label='model backtracked', color='C1', ls='dashed')),
        ]

        fig = plot_twissfigure(model, curves, ylim=(-25, +25))
        fig.savefig(img_base + '-fit-25.png', dpi=400)

        fig = plot_twissfigure(model, curves, ylim=(-40, +40))
        fig.savefig(img_base + '-fit-40.png', dpi=400)

        fig = plot_twissfigure(model, curves, ylim=(-200, +200))
        fig.savefig(img_base + '-fit-200.png', dpi=400)


def backtrack(model, measured, bpms):
    return fit_particle_readouts(model, [
        Readout(
            bpm_name,
            measured[bpm_name]['x'],
            measured[bpm_name]['y'])
        for bpm_name in bpms
        if bpm_name in measured
    ], to='#s')[0][0]


def interpolate(model, measured, bpms):
    init = backtrack(model, measured, bpms)
    model.update_twiss_args(init)
    twiss = model.twiss()
    return twiss.s, twiss.x, twiss.y


def plot_twissfigure(model, curves, ylim=None):
    fig = create_twiss_figure(model, ['sbend'])
    ax_x, ax_y = fig.axes
    ax_x.set_ylabel('x [mm]')
    ax_y.set_ylabel('y [mm]')
    ax_x.axvspan(
        model.elements['gant_rot'].position,
        model.elements['patient'].position,
        color='black', alpha=0.03)
    ax_y.axvspan(
        model.elements['gant_rot'].position,
        model.elements['patient'].position,
        color='black', alpha=0.03)

    for (s, x, y), (args, kwargs) in curves:
        ax_x.plot(s, x * 1000, *args, **kwargs)
        ax_y.plot(s, y * 1000, *args, **kwargs)

    ax_x.legend(loc='lower center', bbox_to_anchor=(0.5, 1.015), ncol=len(curves))
    if ylim:
        clip_ylim(ax_x, *ylim)
        clip_ylim(ax_y, *ylim)
    return fig


def strip_suffix(s, suffix):
    return s[:len(s) - len(suffix)] if s.endswith(suffix) else s


def sort_dict(d, key):
    return {k: d[k] for k in sorted(d, key=key)}


def clip_ylim(ax, lo, hi):
    ax.autoscale()
    a, b = ax.get_ylim()
    ax.set_ylim(max(a, lo), min(b, hi))


def args(*args, **kwargs):
    return args, kwargs


def copy_results():
    shutil.copy(
        '../data/correct/2018-07-03-correct/gantry_p_e1_g0-fit-25.png',
        '../../reports/2019-06-14-madgui/plots/orbit-simple-lim-25.png')
    shutil.copy(
        '../data/correct/2018-07-03-correct/gantry_p_e1_g0-fit-200.png',
        '../../reports/2019-06-14-madgui/plots/orbit-simple-lim-200.png')


def main():
    plot_single_optic_measurements()
    os.makedirs('presentation_plots', exist_ok=True)
    copy_results()


if __name__ == '__main__':
    main()

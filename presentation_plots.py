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
# BPM value evolution

import shutil
from glob import glob

import numpy as np

from madgui.online.orbit import fit_particle_readouts, Readout
from madgui.model.madx import Model
from madgui.util.yaml import load_file
from orm_plot import create_twiss_figure, savefig
from orm_util import parse_errors, filter_errors


# assumed standard error of mean on BPM in addition to statistical error:
BPM_ERR = 0.0003    # [m]


# ====
# DATA
# ====

DATA_PREFIX = '../data/'
HEBT_BPMS = ['h1dg1g', 'h1dg2g', 'h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g',
             'b1dg2g', 'b1dg3g', 't1dg2g', 't1dg1g', 't1df1',
             'b2dg2g', 'b2dg3g', 't2dg2g', 't2dg1g', 't2df1']
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
    ('hht3', 'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/'),
    ('hht3', 'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G4/'),
    ('hht3', 'orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G4-ISO/'),
    ('hht3', 'orm/2019-01-20-orm_measurements/M3-E108-F1-I9-G10-ISO/'),
    ('hht3', 'orm/2019-01-20-orm_measurements/M8-E108-F1-I9-G1-ISO/'),
    ('hht3', 'orm/2019-01-20-quadscan/M8-E108-F1-I9-G1/0-h3qd22/'),
    ('hht3', 'orm/2019-01-20-quadscan/M8-E108-F1-I9-G1/1-b3qd12/'),
    ('hht3', 'orm/2019-01-20-quadscan/M8-E108-F1-I9-G1/2-g3qd22/'),
    ('hht3', 'orm/2019-01-20-quadscan/M8-E108-F1-I9-G1/3-g3qd42/'),
    ('hht3', 'orm/2019-05-11-orm_measurements/M8-E108-F1-I9-G1/'),
    ('hht3', 'orm/2019-05-11-orm_measurements/M8-E108-F1-I9-G1-LIN/'),
    ('hht1', 'orm/2019-06-10-orm_measurements/M6-E108-F1-I9-H1/'),
    ('hht2', 'orm/2019-06-10-orm_measurements/M7-E108-F1-I9-H2/'),
    ('hht3', 'orm/2019-06-10-orm_measurements/M8-E108-F1-I9-G10/'),
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


fit_files = [
#   'fits/a_0_fit_14.yml',
    'fits/a_0_fit_62.yml',
#   'fits/a_2_fit_72.yml',
#   'fits/a_3_fit_63.yml',
#   'fits/b_0_fit.txt',
#   'fits/b_1_fit.txt',
#   'fits/b_2_fit_G1.txt',
#   'fits/b_3_fit_G4.txt',
#   'fits/c_0_fit_41.23.txt',
#   'fits/c_1_fit_15.44.txt',
#   'fits/c_2_fit_40.74.txt',
#   'fits/c_3_fit_24.33.txt',
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
        savefig(fig, img_base + '-fit-25')

        fig = plot_twissfigure(model, curves, ylim=(-40, +40))
        savefig(fig, img_base + '-fit-40')

        fig = plot_twissfigure(model, curves, ylim=(-200, +200))
        savefig(fig, img_base + '-fit-200')


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


def plot_beampos_offset():
    print("\nPlotting beam position / offsets at G3DG5G")
    from orm_util import Analysis
    from backtrack_and_fit_tm import plot_beampos_at
    for seq, folder in orm_measurements_folders:
        record_files = DATA_PREFIX + folder + '*.yml'
        print(record_files)
        ana = Analysis.session(f'../hit_models/{seq}', record_files)
        ana.absolute = False
        ana.measured.stderr = (ana.measured.stderr**2 + BPM_ERR**2)**0.5
        if ana.ensure_monitors_available(['g3dg5g'] + ISOC_BPMS):
            prefix = DATA_PREFIX + folder
            plot_beampos_at(ana, prefix)


def plot_orms():
    print("Plotting ORMs for each session")
    from orm_util import Analysis
    for seq, folder in orm_measurements_folders:
        record_files = DATA_PREFIX + folder + '*.yml'
        name = seq[-2:].upper() + (' (pre gantry)' if seq == 'hht3' else '')
        print(record_files)

        ana = Analysis.session(f'../hit_models/{seq}', record_files)
        ana.absolute = False

        hebt_bpms = [m for m in ana.monitors if m in HEBT_BPMS]
        gant_bpms = [m for m in ana.monitors if m not in HEBT_BPMS
                     and m not in ('t3dg1g', 't3df1')]
        isoc_bpms = [m for m in ana.monitors if m in ISOC_BPMS]

        ana.measured.stderr = (ana.measured.stderr**2 + BPM_ERR**2)**0.5
        ana.model.update_twiss_args(dict(x=0, y=0, px=0, py=0))
        ana.init()
        base_orbits = ana.model_orbits
        ana.info([ana.monitors.index(m) for m in hebt_bpms])
        ana.info([ana.monitors.index(m) for m in gant_bpms])
        if hebt_bpms:
            fig = ana.plot_orm(hebt_bpms)
            suptitle(fig, f"Orbit response: {name}")
            savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-init0-horiz')
        if gant_bpms:
            fig = ana.plot_orm(gant_bpms)
            suptitle(fig, "Orbit response: T3 (gantry)")
            savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-init0-gantry')

        fit_orms = []
        for errfile in fit_files:
            errors = parse_errors(load_file(errfile))
            errors = filter_errors(errors, ana.model)
            ana.apply_errors(errors.keys(), errors.values())
            ana.init()
            fit_orms.append(ana.model_orbits)

        ana.model_orbits = base_orbits
        for ifit, fitted in enumerate(fit_orms):
            if hebt_bpms:
                fig = ana.plot_orm(hebt_bpms, fitted)
                suptitle(fig, f"Orbit response: {name}")
                savefig(fig, DATA_PREFIX + folder[:-1] + f'-orm-f{ifit}-init0-horiz')
            if gant_bpms:
                fig = ana.plot_orm(gant_bpms, fitted)
                suptitle(fig, "Orbit response: T3 (gantry)")
                savefig(fig, DATA_PREFIX + folder[:-1] + f'-orm-f{ifit}-init0-gantry')

        ana = Analysis.session(f'../hit_models/{seq}', record_files)
        ana.absolute = False
        ana.measured.stderr = (ana.measured.stderr**2 + BPM_ERR**2)**0.5
        ana.model.update_twiss_args(dict(x=0, y=0, px=0, py=0))
        ana.init()

        # As it turns out, the ORMs obtained by using different initial
        # conditions are very similar to those assuming x=y=px=py=0. We
        # therefore do currently repeat the fitted error plots here.

        ana.apply_errors([], [])
        if len(hebt_bpms) >= 2:
            ana.backtrack(hebt_bpms)
            if hebt_bpms:
                fig = ana.plot_orm(hebt_bpms)
                suptitle(fig, f"Orbit response: {name}")
                savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-backfit-horiz')
            if gant_bpms:
                fig = ana.plot_orm(gant_bpms)
                suptitle(fig, "Orbit response: T3 (gantry)")
                savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-backfit-gantry')

        if len(isoc_bpms) >= 2:
            ana.backtrack(isoc_bpms)
            if hebt_bpms:
                fig = ana.plot_orm(hebt_bpms)
                suptitle(fig, f"Orbit response: {name}")
                savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-backisoc-horiz')
            if gant_bpms:
                fig = ana.plot_orm(gant_bpms)
                suptitle(fig, "Orbit response: T3 (gantry)")
                savefig(fig, DATA_PREFIX + folder[:-1] + '-orm-backisoc-gantry')


def suptitle(fig, title):
    fig.suptitle(title, x=0.04, y=0.98,
                 horizontalalignment='left',
                 verticalalignment='top')


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
    dest = '../../reports/2019-06-14-madgui/fig/'
    shutil.copy(
        '../data/correct/2018-07-03-correct/gantry_p_e1_g0-fit-25.pdf',
        dest + 'orbit-simple-lim-25.pdf')
    shutil.copy(
        '../data/correct/2018-07-03-correct/gantry_p_e1_g0-fit-200.pdf',
        dest + 'orbit-simple-lim-200.pdf')
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/beampos-g3dg5g.pdf',
        dest)
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/beampos_vs_optic.pdf',
        dest)
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1-orm-backfit-horiz.pdf',
        dest + 'orm-backfit-horiz.pdf')
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1-orm-backfit-gantry.pdf',
        dest + 'orm-backfit-gantry.pdf')
    # fitted
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1-orm-f0-init0-horiz.pdf',
        dest + 'orm-fitted-horiz.pdf')
    shutil.copy(
        '../data/orm/2018-10-20-orm_measurements/M8-E108-F1-I9-G1-orm-f0-init0-gantry.pdf',
        dest + 'orm-fitted-gantry.pdf')
    # new results
    shutil.copy(
        '../data/orm/2019-06-10-orm_measurements/M6-E108-F1-I9-H1-orm-init0-horiz.pdf',
        dest + 'orm-h1.pdf')
    shutil.copy(
        '../data/orm/2019-06-10-orm_measurements/M7-E108-F1-I9-H2-orm-init0-horiz.pdf',
        dest + 'orm-h2.pdf')
    shutil.copy(
        '../data/orm/2019-06-10-orm_measurements/M8-E108-F1-I9-G10-orm-init0-horiz.pdf',
        dest + 'orm-h3-horiz.pdf')
    shutil.copy(
        '../data/orm/2019-06-10-orm_measurements/M8-E108-F1-I9-G10-orm-init0-gantry.pdf',
        dest + 'orm-h3-gantry.pdf')


def main():
    from madgui.core.app import init_app
    init_app(['madgui'])
    plot_single_optic_measurements()
    plot_beampos_offset()
    plot_orms()
    copy_results()


if __name__ == '__main__':
    main()

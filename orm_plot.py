"""
Plot utility functions for showing results of the ORM analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from madgui.plot.twissfigure import (
    plot_element_indicators, with_outline, ELEM_STYLES)


def plot_orm(model, measured, orbits, monitors):
    fig = plt.figure(1)
    fig.clf()

    monitors = [m for m in monitors if m in measured.monitors]

    model_orm = orbits[:, :, 1:] - orbits[:, :, [0]]
    measured_orm = measured.orbits[:, :, 1:] - measured.orbits[:, :, [0]]
    measured_errors = (measured.stderr[:, :, 1:] ** 2 +
                       measured.stderr[:, :, [0]] ** 2) ** 0.5

    for iy, name in enumerate("xy"):
        ax = fig.add_subplot(2, 1, 1+iy)
        ax.set_ylabel(f'$\Delta {name}$ [mm]')

        num_entries = []
        ydata_model = []
        ydata_measured = []
        errors = []

        for mon in monitors:
            idx = measured.monitors.index(mon)
            entry_model = model_orm[idx, iy] * 1000
            entry_measured = measured_orm[idx, iy] * 1000
            entry_errors = measured_errors[idx, iy] * 1000
            defined_region = ~np.isnan(entry_measured)
            ydata_model.extend(entry_model[defined_region])
            ydata_measured.extend(entry_measured[defined_region])
            errors.extend(entry_errors[defined_region])
            num_entries.append(round(sum(defined_region)))

        xdata = range(len(ydata_model))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.plot(xdata, ydata_model, '-', color='C0', label="model")
        ax.errorbar(xdata, ydata_measured, errors, fmt='.',
                    color='C1', label="measured")

        cumsum = np.cumsum(num_entries)
        for pos in cumsum[:-1]:
            ax.axvline(pos + 0.5, linestyle='--', color='k', linewidth=0.5)

    for mon, pos in zip(monitors, cumsum):
        ax.text(pos-0.5, 0.05, mon,
                horizontalalignment='right',
                verticalalignment='bottom',
                rotation=90,
                alpha=0.6,
                transform=ax.get_xaxis_transform())

    fig.axes[-1].set_xlabel("entry")
    fig.axes[0].legend(
        loc='lower right', bbox_to_anchor=(1.0, 1.05), ncol=2)
    fig.suptitle("Orbit response")
    return fig


def make_orbit_plots(
        model, measured, orbits, optics,
        save_to=None, base_orbit=None):

    for i, (optic, tw) in enumerate(zip(optics, orbits)):
        fig = create_twiss_figure(model)
        plot_orbit(fig, model, i, tw, measured,
                   base_orbit=base_orbit and base_orbit[i])
        fig.suptitle(f"Orbit for optic #{i}")
        knob = next(iter(optic), (None, None))[0]
        savefig(fig, save_to and f'{save_to}-orbit-{i}-{knob}', save_to)
        plt.clf()


def make_monitor_plots(
        monitor_subset, model, measured, model_orbits, comment="Response",
        save_to=None, base_orm=None):
    for index, monitor in enumerate(measured.monitors):
        if monitor in monitor_subset:
            fig = create_twiss_figure(model)
            plot_monitor_response(
                fig, monitor,
                model, measured, base_orm, model_orbits, comment)
            fig.suptitle("Orbit response at {monitor}")
            savefig(fig, save_to and f'{save_to}-mon-{index}-{monitor}')
            plt.clf()


def make_steerer_plots(
        steerer_subset, model, measured, model_orbits, comment="Response",
        save_to=None, base_orm=None):
    for index, steerer in enumerate(measured.steerers):
        if steerer in steerer_subset:
            fig = create_twiss_figure(model)
            plot_steerer_response(
                fig, steerer,
                model, measured, base_orm, model_orbits, comment)
            fig.suptitle("Orbit response due to {steerer}")
            savefig(fig, save_to and f'{save_to}-ste-{index}-{steerer}')
            plt.clf()


def plot_monitor_response(
        fig, monitor, model, measured, base_orm, model_orbits, comment):
    xpos = [model.elements[elem].position for elem in measured.steerers]
    i = measured.monitors.index(monitor)
    lines = []

    orbits = response_matrix(measured.orbits)
    model_orbits = response_matrix(model_orbits)
    base_orm = response_matrix(base_orm)
    stderr = measured.stderr
    if stderr is not None:
        stderr = (stderr[:, :, 1:]**2 + stderr[:, :, [0]]**2)**0.5

    for j, ax, axes in zip(range(2), "xy", fig.axes):

        axes.set_title(ax)
        axes.set_xlabel(r"steerer position [m]")
        if ax == 'x':
            axes.set_ylabel(r"orbit response $\Delta x/\Delta \phi$ [mm/mrad]")
        else:
            axes.yaxis.tick_right()

        axes.errorbar(
            xpos,
            orbits[i, j, :].flatten(),
            stderr[i, j, :].flatten(),
            label=ax + " measured")

        lines.append(axes.plot(
            xpos,
            model_orbits[i, j, :].flatten(),
            label=ax + " model"))

        if base_orm is not None:
            axes.plot(
                xpos,
                base_orm[i, j, :].flatten(),
                label=ax + " base model")

        axes.legend()

    fig.suptitle("{1}: {0}".format(monitor, comment))
    return lines


def plot_steerer_response(
        fig, steerer, model, measured, base_orm, model_orbits, comment):
    xpos = [model.elements[elem].position for elem in measured.monitors]
    i = measured.steerers.index(steerer)
    lines = []

    orbits = response_matrix(measured.orbits)
    model_orbits = response_matrix(model_orbits)
    base_orm = response_matrix(base_orm)
    stderr = measured.stderr
    if stderr is not None:
        stderr = (stderr[:, :, 1:]**2 + stderr[:, :, [0]]**2)**0.5

    for j, ax, axes in zip(range(2), "xy", fig.axes):
        axes.set_title(ax)
        axes.set_xlabel(r"Position [m]")
        if ax == 'x':
            axes.set_ylabel(r"orbit response $\Delta x/\Delta \phi$ [mm/mrad]")
        else:
            axes.yaxis.tick_right()

        axes.errorbar(
            xpos,
            orbits[:, j, i].flatten(),
            stderr[:, j, i].flatten(),
            label=ax + " measured")

        lines.append(axes.plot(
            xpos,
            model_orbits[:, j, i].flatten(),
            label=ax + " model"))

        if base_orm is not None:
            axes.plot(
                xpos,
                base_orm[:, j, i].flatten(),
                label=ax + " base model")

        axes.legend()

    fig.suptitle("{1}: {0}".format(steerer, comment))
    return lines


def response_matrix(orbits):
    return None if orbits is None else orbits[:, :, 1:] - orbits[:, :, [0]]


def plot_orbit(fig, model, i, twiss, measured, base_orbit):

    xpos = [model.elements[elem].position for elem in measured.monitors]
    orbit = measured.orbits[:, :, i]
    error = measured.stderr[:, :, i]

    for j, ax, axes in zip(range(2), "xy", fig.axes):

        axes.set_ylabel(rf"${ax}$ [mm]")

        axes.errorbar(xpos, orbit[:, j] * 1000, error[:, j] * 1000,
                      fmt='o', label="measured", markersize=5)
        axes.plot(twiss.s, twiss[ax]*1000, label="model",
                  **with_outline({}))
        if base_orbit is not None:
            axes.plot(
                base_orbit.s, base_orbit[ax], label=ax + " base_orbit")

        if ax == 'x':
            axes.legend(loc='upper right')

    fig.suptitle("orbit")


def create_twiss_figure(model, elem_types=None):
    fig = plt.figure(1)
    fig.clf()
    elem_styles = ELEM_STYLES
    if elem_types is not None:
        elem_styles = {k: v for k, v in elem_styles.items() if k in elem_types}
    axx = fig.add_subplot(2, 1, 1)
    axy = fig.add_subplot(2, 1, 2, sharex=axx)
    axy.set_xlabel("Position $s$ [m]")
    for ax, axis in zip((axx, axy), "xy"):
        ax.x_name = ['s']
        ax.y_name = [axis]
        plot_element_indicators(
            ax, model.elements, elem_styles, effects=background_style)
        ax.grid(True, axis='y')
    return fig


def background_style(style):
    alpha = 1 if style.get('alpha') == 1 else 0.25
    return dict(style, alpha=alpha)


def savefig(fig, to):
    if to is None:
        plt.show(fig)
    else:
        fig.savefig(to + '.png', dpi=400, bbox_inches='tight')
        fig.savefig(to + '.pdf')

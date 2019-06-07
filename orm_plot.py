"""
Plot utility functions for showing results of the ORM analysis.
"""

import matplotlib.pyplot as plt
from madgui.plot.twissfigure import plot_element_indicators, with_outline


def make_orbit_plots(
        model, measured, orbits, optics,
        save_to=None, base_orbit=None):

    for i, (optic, tw) in enumerate(zip(optics, orbits)):
        fig = create_twiss_figure(model)
        plot_orbit(fig, model, i, tw, measured,
                   base_orbit=base_orbit and base_orbit[i])
        fig.suptitle(f"Orbit for optic #{i}")
        if save_to is None:
            plt.show()
        else:
            knob = next(iter(optic), (None, None))[0]
            plt.savefig(f'{save_to}-orbit-{i}-{knob}.png', dpi=300)
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
            if save_to is None:
                plt.show()
            else:
                plt.savefig(f'{save_to}-mon-{index}-{monitor}.png', dpi=300)
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
            if save_to is None:
                plt.show()
            else:
                plt.savefig(f'{save_to}-ste-{index}-{steerer}.png', dpi=300)
            plt.clf()


def plot_monitor_response(
        fig, monitor, model, measured, base_orm, model_orbits, comment):
    xpos = [model.elements[elem].position for elem in measured.steerers]
    i = measured.monitors.index(monitor)
    lines = []

    orbits = response_matrix(measured.orbits)
    model_orbits = response_matrix(model_orbits)
    base_orm = response_matrix(base_orm)
    stddev = measured.stddev
    if stddev is not None:
        stddev = (stddev[:, :, 1:]**2 + stddev[:, :, [0]]**2)**0.5

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
            stddev[i, j, :].flatten(),
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
    stddev = measured.stddev
    if stddev is not None:
        stddev = (stddev[:, :, 1:]**2 + stddev[:, :, [0]]**2)**0.5

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
            stddev[:, j, i].flatten(),
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
    error = measured.stddev[:, :, i]

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


def create_twiss_figure(model):
    fig = plt.figure(1)
    axx = fig.add_subplot(2, 1, 1)
    axy = fig.add_subplot(2, 1, 2, sharex=axx)
    axy.set_xlabel("Position $s$ [m]")
    for ax, axis in zip((axx, axy), "xy"):
        ax.x_name = ['s']
        ax.y_name = [axis]
        plot_element_indicators(ax, model.elements, effects=background_style)
        ax.grid(True, axis='y')
    return fig


def background_style(style):
    return dict(style, alpha=0.2)

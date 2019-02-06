"""
Plot utility functions for showing results of the ORM analysis.
"""

import matplotlib.pyplot as plt


def make_orbit_plots(
        model, measured, orbits, optics,
        save_to=None, base_orbit=None):

    for i, (optic, tw) in enumerate(zip(optics, orbits)):
        fig = plt.figure(1)
        plot_orbit(fig, model, i, tw, measured,
                   base_orbit=base_orbit and base_orbit[i])
        if save_to is None:
            plt.show()
        else:
            knob = next(iter(optic), (None, None))[0]
            plt.savefig('{}-orbit-{}-{}.png'.format(save_to, i, knob))
        plt.clf()


def make_monitor_plots(
        monitor_subset, model, measured, model_orbits, comment="Response",
        save_to=None, base_orm=None):
    for index, monitor in enumerate(measured.monitors):
        if monitor in monitor_subset:
            plot_monitor_response(
                plt.figure(1), monitor,
                model, measured, base_orm, model_orbits, comment)
            if save_to is None:
                plt.show()
            else:
                plt.savefig('{}-mon-{}-{}.png'.format(save_to, index, monitor))
            plt.clf()


def make_steerer_plots(
        steerer_subset, model, measured, model_orbits, comment="Response",
        save_to=None, base_orm=None):
    for index, steerer in enumerate(measured.steerers):
        if steerer in steerer_subset:
            plot_steerer_response(
                plt.figure(1), steerer,
                model, measured, base_orm, model_orbits, comment)
            if save_to is None:
                plt.show()
            else:
                plt.savefig('{}-ste-{}-{}.png'.format(save_to, index, steerer))
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

    for j, ax in enumerate("xy"):
        axes = fig.add_subplot(1, 2, 1+j)
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

    for j, ax in enumerate("xy"):
        axes = fig.add_subplot(1, 2, 1+j)
        axes.set_title(ax)
        axes.set_xlabel(r"monitor position [m]")
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

    for j, ax in enumerate("xy"):
        axes = fig.add_subplot(1, 2, 1+j)
        axes.set_title(ax)
        axes.set_xlabel(r"monitor position [m]")
        if ax == 'x':
            axes.set_ylabel(r"orbit response $\Delta x/\Delta \phi$ [mm/mrad]")
        else:
            axes.yaxis.tick_right()

        axes.errorbar(xpos, orbit[:, j], error[:, j],
                      fmt='v-', label=ax + " measured")
        axes.plot(twiss.s, twiss[ax], label=ax + " model")
        if base_orbit is not None:
            axes.plot(base_orbit.s, base_orbit[ax], label=ax + " base_orbit")
        axes.legend()

    fig.suptitle("orbit")

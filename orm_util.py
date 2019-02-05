from contextlib import contextmanager

import matplotlib.pyplot as plt
import numpy as np

from cpymad.madx import TwissFailed

import madgui.util.yaml as yaml
from madgui.util.fit import reduced_chisq, fit
from madgui.online.orbit import fit_particle_readouts, Readout
from madgui.model.errors import Ealign, Efcomp, apply_errors


class OrbitResponse:

    def __init__(self, strengths, records, steerers, monitors, optics):
        self.monitors = monitors
        self.steerers = steerers
        self.optics = optics
        self.strengths = strengths
        orbits = {
            key: np.mean(shots, axis=0)
            for key, shots in records.items()
        }
        errors = {
            key: np.var(shots, axis=0, ddof=1)
            for key, shots in records.items()
        }
        counts = {
            key: len(shots)
            for key, shots in records.items()
        }

        nans = [np.nan, np.nan]
        # TODO: we should use the "base optic" before the "steerer"
        self.orbits = np.dstack([
            np.vstack([
                orbits.get((monitor, optic), nans)
                for monitor in monitors
            ])
            for optic in optics
        ])
        self.variance = np.dstack([
            np.vstack([
                errors.get((monitor, optic), nans)
                for monitor in monitors
            ])
            for optic in optics
        ])
        self.stddev = self.variance ** 0.5 / np.dstack([
            np.vstack([
                counts.get((monitor, optic), 1) ** 0.5
                for monitor in monitors
            ])
            for optic in optics
        ])
        self.counts = counts

    @classmethod
    def load(cls, model, filenames):
        strengths = {}
        records = {}
        for s, r in map(load_record_file, filenames):
            # TODO: should handle different base_optics in different files…
            strengths.update(s)
            for k, v in r.items():
                records.setdefault(k, []).extend(v)
        monitors = {mon.lower() for mon, _ in records}
        optics = {optic for _, optic in records}
        knob_elems = map_knobs_to_elements(model)
        steerers = {knob_elems[next(iter(optic))[0]].name
                    for optic in optics if optic}
        return cls(strengths, records, steerers,
                   sorted(monitors, key=model.elements.index),
                   sorted_optics(model, optics))


def map_knobs_to_elements(model):
    """Build a dictionary that maps knob name (lowercase) to element."""
    knob_elems = {}
    for elem in model.elements:
        for knob in model.get_elem_knobs(elem):
            knob_elems.setdefault(knob.lower(), elem)
    return knob_elems


def sorted_optics(model, optics):
    """Sort lists of knob/value pairs by element index and knob value."""
    knob_elems = map_knobs_to_elements(model)
    elem_index = {knob: elem.index for knob, elem in knob_elems.items()}
    def item_key(item):
        knob, value = item
        return (elem_index.get(knob, -1), knob, value)
    optics = [tuple(sorted(optic, key=item_key)) for optic in optics]
    return sorted(optics, key=lambda optic: [
        item_key(item) for item in optic])


def load_record_file(filename):
    """Load a YAML file in the format produced by the madgui dialog at
    Online control -> ORM measurement."""
    data = yaml.load_file(filename)
    strengths = data['model']
    records = {
        (monitor, tuple(record['optics'].items())): [
            shot[monitor][:2]
            for shot in record['shots']
        ]
        for record in data['records']
        for monitor in data['monitors']
    }
    return strengths, records


def fit_init_orbit(model, measured, fit_monitors):
    (twiss_init, chisq, singular), curve = fit_particle_readouts(model, [
        Readout(monitor, *measured.orbits[index, :, 0])
        for monitor in fit_monitors
        for index in [measured.monitors.index(monitor.lower())]
    ])
    return twiss_init


class Analysis:

    def __init__(self, model, measured):
        self.model = model
        self.measured = measured
        self.monitors = measured.monitors
        self.steerers = measured.steerers
        self.optics = measured.optics
        self.errors = []
        self.values = []

    def init(self, strengths=None):
        print("INITIAL")
        if strengths is None:
            strengths = self.measured.strengths
        self.model.update_globals(strengths.items())
        self.model_orbits = self.compute_model_orbits()
        sel = self.get_selected_monitors(self.monitors)
        self.info(sel)

    def info(self, sel=None, ddof=0):
        if sel is None:
            sel = slice(None)
        measured = self.measured
        model_orbits = self.model_orbits
        stddev = measured.stddev
        print("red χ² =", reduced_chisq(
            ((measured.orbits - model_orbits) / stddev)[sel], ddof))
        print("    |x =", reduced_chisq(
            ((measured.orbits - model_orbits) / stddev)[sel][:, 0, :], ddof))
        print("    |y =", reduced_chisq(
            ((measured.orbits - model_orbits) / stddev)[sel][:, 1, :], ddof))

    def apply_errors(self, errors, values):
        self.errors[:0] = errors
        self.values[:0] = values

    def compute_model_orbits(self, errors=(), values=()):
        errs = list(errors) + self.errors
        vals = list(values) + self.values
        idx = [self.model.elements.index(m) for m in self.monitors]
        return np.dstack([
            self._get_orbit(optic, errs, vals)
            .dframe(('x', 'y'))
            for optic in self.optics
        ])[idx]

    def _get_orbit(self, optic, errs, vals):
        return get_orbit(self.model, optic, errs, vals)

    def get_selected_monitors(self, selected):
        return [self.monitors.index(m.lower()) for m in selected]

    def plot_monitors(self, select=None, save_to=None, base_orm=None):
        if select is None:
            select = self.monitors
        print("plotting monitors: {}".format(" ".join(select)))
        make_monitor_plots(
            select, self.model, self.measured, self.model_orbits,
            save_to=save_to, base_orm=base_orm)

    def plot_steerers(self, select=None, save_to=None, base_orm=None):
        if select is None:
            select = self.steerers
        print("plotting steerers: {}".format(" ".join(select)))
        make_steerer_plots(
            select, self.model, self.measured, self.model_orbits,
            save_to=save_to, base_orm=base_orm)

    def plot_orbit(self, save_to=None, base_orbit=None):

        optics = self.optics
        orbits = [
            self._get_orbit(optic, self.errors, self.values)
            .dframe(('s', 'x', 'y'))
            for optic in optics
        ]

        for i, (optic, tw) in enumerate(zip(optics, orbits)):
            fig = plt.figure(1)
            plot_orbit(fig, self.model, i, tw, self.measured,
                       base_orbit=base_orbit and base_orbit[i])
            if save_to is None:
                plt.show()
            else:
                knob = next(iter(optics), (None, None))[0]
                plt.savefig('{}-orbit-{}-{}.png'.format(save_to, i, knob))
            plt.clf()

        return orbits

    def backtrack(self, monitors):
        print("TWISS INIT")
        twiss_args = fit_init_orbit(self.model, self.measured, monitors)
        self.model.update_twiss_args(twiss_args)
        self.model_orbits = self.compute_model_orbits()
        return twiss_args

    def fit(self, errors, monitors, delta=1e-4,
            mode='xy', iterations=50, bounds=None, fourier=False,
            tol=1e-8, use_stddev=True, save_to=None, **kwargs):

        model = self.model
        measured = self.measured
        stddev = measured.stddev if use_stddev else 1
        err_names = ', '.join(map(repr, errors))

        print("====================")
        print("FIT:", ', '.join(monitors or self.monitors))
        print("VIA:", err_names)

        sel = self.get_selected_monitors(monitors or self.monitors)
        inv = sorted(set(range(len(self.monitors))) - set(sel))
        model.madx.eoption(add=True)

        def callback(state):
            print("")
            print("----------------------")
            print("nit    =", state.nit)
            print("Errors :", err_names)
            print("ΔX     =", state.dx)
            print("X_tot  =", state.x)
            print(":: (fit) ::")
            self.info(sel)
            if inv:
                print(":: (elsewhere) ::")
                self.info(inv)
                print(":: (overall) ::")
                self.info()
            print("----------------------")

        dims = [i for i, c in enumerate("xy") if c in mode]

        def objective(values):
            try:
                print(".", end='', flush=True)
                self.model_orbits = self.compute_model_orbits(errors, values)
            except TwissFailed:
                return np.array([1e8])
            obj = ((self.model_orbits - measured.orbits) / stddev)[sel][:, dims, :]
            if fourier:
                obj = np.fft.rfft(obj, axis=0)
                obj = np.array([
                    np.real(obj),
                    np.imag(obj),
                ]).transpose((1, 2, 3, 0))
            return obj

        x0 = np.zeros(len(errors))
        result = fit(
            objective, x0, tol=tol,
            delta=delta, iterations=iterations, callback=callback, **kwargs)
        print(result.message)
        self.apply_errors(errors, result.x)

        if save_to is not None:
            text = '\n'.join(
                '{!r}: {}'.format(err, val)
                for err, val in zip(errors, result.x))
            with open(save_to, 'wt') as f:
                f.write(text)

        return result

    @classmethod
    @contextmanager
    def app(cls, model_file, record_files):
        from madgui.core.app import init_app
        from madgui.core.session import Session
        from madgui.core.config import load as load_config
        from glob import glob

        init_app(['madgui'])

        if isinstance(record_files, str):
            record_files = glob(record_files)

        config = load_config(isolated=True)
        with Session(config) as session:
            session.load_model(
                model_file,
                stdout=False)
            model = session.model()
            measured = OrbitResponse.load(model, record_files)
            yield cls(model, measured)


def get_orbit(model, optic, errors, values, **twiss_args):
    """Get x, y vectors, with specified errors."""
    twiss_args.setdefault('table', 'orm_tmp')
    madx = model.madx
    madx.command.select(flag='interpolate', clear=True)
    with model.undo_stack.rollback():
        model.update_globals(optic)
        with apply_errors(model, errors, values):
            return madx.twiss(**model._get_twiss_args(**twiss_args))


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


ERR_ATTR = {
    'sbend': ['angle', 'e1', 'e2', 'k0', 'hgap', 'fint'],
    'quadrupole': ['k1', 'k1s'],
    'hkicker': ['kick', 'tilt'],
    'vkicker': ['kick', 'tilt'],
    'srotation': ['angle'],
}

ERR_EALIGN = ['dx', 'dy', 'ds', 'dpsi', 'dphi', 'dtheta']


# scaling: monitor
# scaling: kicker

def get_elem_ealign(model, name, attrs=ERR_EALIGN, delta=1e-3):
    return [
        Ealign({'range': name}, attr, delta)
        for attr in attrs
    ]


def get_elem_efcomp(model, name, delta=1e-3):
    elem = model.elements[name]
    kwargs = dict(order=None, radius=None)
    if elem.base_name == 'sbend':
        return [Efcomp({'range': name}, 'dkn', [delta], **kwargs)]
    if elem.base_name == 'quadrupole':
        return [
            Efcomp({'range': name}, 'dkn', [0, delta], **kwargs),
            Efcomp({'range': name}, 'dks', [0, delta], **kwargs),
        ]
    return []

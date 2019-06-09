import itertools

import numpy as np

from cpymad.madx import TwissFailed

import madgui.util.yaml as yaml
from madgui.util.fit import reduced_chisq, fit
from madgui.online.orbit import fit_particle_readouts, Readout
from madgui.model.errors import apply_errors, parse_error
from madgui.util.undo import UndoStack


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
        # position measurements, stacked as [monitor × axis × optic]
        self.orbits = np.dstack([
            np.vstack([
                orbits.get((monitor, optic), nans)
                for monitor in monitors
            ])
            for optic in optics
        ])
        # variance of the distribution: σ² = Σ(x - mu_x)² / (N-1)
        # (using bessel's correction N-1)
        self.variance = np.dstack([
            np.vstack([
                errors.get((monitor, optic), nans)
                for monitor in monitors
            ])
            for optic in optics
        ])
        # standard error of mean: SEM = σ/√N
        self.stderr = self.variance ** 0.5 / np.dstack([
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
        return cls(strengths, records,
                   sorted(steerers, key=model.elements.index),
                   sorted(monitors, key=model.elements.index),
                   sorted_optics(model, optics))

    def filter_optics(self, indices):
        self.optics = [self.optics[i] for i in indices]
        self.orbits = self.orbits[:, :, indices]
        self.variance = self.variance[:, :, indices]
        self.stderr = self.stderr[:, :, indices]

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
    """
    Load a YAML file in the format produced by one of the madgui dialogs:

        - Online control -> ORM measurement
        - Online control -> Beam diagnostic -> Offsets -> Calibrate
    """
    data = yaml.load_file(filename)
    if 'base_optics' in data:
        # exported by `Online -> Beam diagnostic -> Offsets -> Calibrate`…
        data = _convert_orm_export(data)
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


def _convert_orm_export(data):
    """
    Convert data as exported by the madgui dialogs at:

    from:       Online control -> Beam diagnostic -> Offsets -> Calibrate
    to:         Online control -> ORM measurement
    """
    return {
        'model': data['base_optics'],
        'monitors': data['monitors'],
        'knobs': {knob for optic in data['optics'] for knob in optic},
        'sequence': 'hht3',
        'records': [
            {
                'optics': optic,
                'shots': [
                    {monitor: [r['posx'], r['posy'], r['envx'], r['envy']]
                     for monitor in data['monitors']
                     for r in [shot['readout'][monitor]]}
                    for shot in group
                ],
            }
            for optic, group in itertools.groupby(
                    data['records'], key=lambda r: r['optics'])
        ],
    }


class Analysis:

    def __init__(self, model, measured):
        self.model = model
        self.measured = measured
        self.monitors = measured.monitors
        self.steerers = measured.steerers
        self.optics = measured.optics
        self.errors = []
        self.values = []
        self._init_twiss = {}
        self.absolute = True

    def ensure_monitors_available(self, monitors):
        """Pick optics for which we have measured all of the given BPMs."""
        if not all(m in self.monitors for m in monitors):
            return False
        involved_bpms = [self.monitors.index(m) for m in monitors]
        usable_optics = [
            i_optic for i_optic in range(len(self.optics))
            if not np.any(np.isnan(
                self.measured.orbits[involved_bpms][:, :, i_optic]))
        ]
        self.measured.filter_optics(usable_optics)
        self.optics = self.measured.optics
        return bool(self.optics)

    def setup_backtracking(self, final_orbits):
        """Reverse sequence, set the initial coordinates for each optic from
        the given final coordinates."""
        self._init_twiss = {
            optic: {'x': -orbit['x'], 'px': orbit['px'],
                    'y': orbit['y'], 'py': -orbit['py']}
            for optic, orbit in final_orbits.items()
        }
        self.model.reverse()
        self.model.update_twiss_args(self._init_twiss.get((), {}))
        self.measured.orbits[:, 0, :] *= -1
        #self.measured.stderr[:, :, :] = 1e-4

    def extrapolate(self, monitors, to='#e'):
        """Extrapolate x/px/y/py for all known optics from the measurements
        at the given monitors."""
        self.ensure_monitors_available(monitors)
        return {
            optic: extrapolate_orbit(
                self.measured, i, optic, self.model, monitors, to=to)
            for i, optic in enumerate(self.optics)
        }

    def init(self, strengths=None):
        print("INITIAL")
        if strengths is None:
            strengths = self.measured.strengths
        self.model.update_globals(strengths.items())
        self.model_orbits = self.compute_model_orbits()
        sel = self.get_selected_monitors(self.monitors)
        self.info(sel)

    def objective(self, model, use_stderr=True):
        measured = self.measured.orbits
        stderr = self.measured.stderr
        if not self.absolute:
            model = model[:, :, 1:] - model[:, :, [0]]
            measured = measured[:, :, 1:] - measured[:, :, [0]]
            stderr = (stderr[:, :, 1:]**2 + stderr[:, :, [0]]**2) ** 0.5
        return (measured - model) / (stderr if use_stderr else 1)

    def info(self, sel=None, ddof=0):
        if sel is None:
            sel = slice(None)
        objective = self.objective(self.model_orbits)
        print("red χ² =", reduced_chisq(objective[sel], ddof))
        print("    |x =", reduced_chisq(objective[sel][:, 0, :], ddof))
        print("    |y =", reduced_chisq(objective[sel][:, 1, :], ddof))

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
        return get_orbit(
            self.model, optic, errs, vals,
            **self._init_twiss.get(optic, {}))

    def get_selected_monitors(self, selected):
        return [self.monitors.index(m.lower()) for m in selected]

    def plot_orm(self, monitors=None):
        if monitors is None:
            monitors = self.monitors
        from orm_plot import plot_orm
        return plot_orm(self.model, self.measured, self.model_orbits, monitors)

    def plot_monitors(self, select=None, save_to=None, base_orm=None):
        if select is None:
            select = self.monitors
        print("plotting monitors: {}".format(" ".join(select)))
        from orm_plot import make_monitor_plots
        make_monitor_plots(
            select, self.model, self.measured, self.model_orbits,
            save_to=save_to, base_orm=base_orm)

    def plot_steerers(self, select=None, save_to=None, base_orm=None):
        if select is None:
            select = self.steerers
        print("plotting steerers: {}".format(" ".join(select)))
        from orm_plot import make_steerer_plots
        make_steerer_plots(
            select, self.model, self.measured, self.model_orbits,
            save_to=save_to, base_orm=base_orm)

    def plot_orbit(self, save_to=None, base_orbit=None):
        orbits = [
            self._get_orbit(optic, self.errors, self.values)
            .dframe(('s', 'x', 'y'))
            for optic in self.optics
        ]
        from orm_plot import make_orbit_plots
        make_orbit_plots(
            self.model, self.measured, orbits, self.optics,
            save_to=save_to, base_orbit=base_orbit)
        return orbits

    def backtrack(self, monitors):
        print("TWISS INIT")
        twiss_args = extrapolate_orbit(
            self.measured, 0, self.optics[0], self.model, monitors, '#s')
        self.model.update_twiss_args(twiss_args)
        self.model_orbits = self.compute_model_orbits()
        return twiss_args

    def fit(self, errors, monitors=None, delta=1e-4,
            mode='xy', iterations=50, bounds=None, fourier=False,
            tol=1e-8, use_stderr=True, save_to=None, **kwargs):

        if isinstance(errors, dict):
            x0 = np.array(list(errors.values()))
            errors = list(errors.keys())
        else:
            x0 = np.zeros(len(errors))

        model = self.model
        err_names = ', '.join(map(repr, errors))

        if monitors is None:
            monitors = self.monitors

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
            obj = self.objective(self.model_orbits, use_stderr)[sel][:, dims, :]
            if fourier:
                obj = np.fft.rfft(obj, axis=0)
                obj = np.array([
                    np.real(obj),
                    np.imag(obj),
                ]).transpose((1, 2, 3, 0))
            return obj

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
    def app(cls, model_file, record_files):
        from madgui.core.app import init_app
        init_app(['madgui'])
        return cls.session(model_file, record_files)

    @classmethod
    def session(cls, model_file, record_files):
        from madgui.core.session import Session
        from madgui.core.config import load as load_config
        from glob import glob
        if isinstance(record_files, str):
            record_files = glob(record_files, recursive=True)

        config = load_config(isolated=True)
        session = Session(config)
        session.load_model(model_file, stdout=False, undo_stack=UndoStack())
        model = session.model()
        measured = OrbitResponse.load(model, record_files)
        analysis = cls(model, measured)
        analysis.session = session
        return analysis


def get_orbit(model, optic, errors, values, **twiss_args):
    """Get x, y vectors, with specified errors."""
    twiss_args.setdefault('table', 'orm_tmp')
    madx = model.madx
    madx.command.select(flag='interpolate', clear=True)
    with model.undo_stack.rollback():
        model.update_globals(optic)
        with apply_errors(model, errors, values):
            return madx.twiss(**model._get_twiss_args(**twiss_args))


def extrapolate_orbit(measured, i_optic, optic, model, from_monitors, to='#e'):
    """Extrapolate particle position/momentum from the position measurements
    of the given BPMs ``from_monitors``."""
    # TODO: in the more general case, we would also need to set the strengths
    # corresponding to i_optic
    with model.undo_stack.rollback():
        model.update_globals(optic)
        return fit_particle_readouts(model, [
            Readout(monitor, *measured.orbits[index, :, i_optic])
            for monitor in from_monitors
            for index in [measured.monitors.index(monitor.lower())]
        ], to=to)[0][0]


def parse_errors(names):
    if isinstance(names, dict):
        return {parse_error(k): v for k, v in names.items()}
    else:
        return [parse_error(k) for k in names]

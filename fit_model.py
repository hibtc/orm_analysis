import sys
from madgui.model.errors import parse_error
from madgui.model.orm import Analysis
from scipy.optimize import Bounds


def parse_errors(names):
    return list(map(parse_error, names))


record_files = (
    sys.argv[1:] or
    '../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1/*.yml')

with Analysis.app('../hit_models/hht3', record_files) as ana:

    ana.init()
#   ana.plot_monitors(save_to='plots/0-initial')
#   ana.plot_steerers(save_to='plots/0-initial')

    #ana.backtrack(['h1dg1g', 'h1dg2g', 'h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g'])
    ana.backtrack(['h2dg2g', 'h3dg3g', 'b3dg2g', 'b3dg3g'])
    #ana.backtrack(['g3dg3g', 'g3dg5g'])
#   ana.plot_monitors(save_to='plots/1-backtrack')
#   ana.plot_steerers(save_to='plots/1-backtrack')

    init_tw = parse_errors(['x', 'y', 'px', 'py'])
    kl_first = parse_errors([
        'δkl_g3qd11',
        'δkl_g3qd12',
        'δkl_g3qd21',
        'δkl_g3qd22',
    ])
    kl_second = parse_errors([
        'δkl_g3qd31',
        'δkl_g3qd32',
        'δkl_efg_g3qd41',
        'δkl_efg_g3qd42',
    ])
    dx_first = parse_errors([
        'g3qd11<dx>',
        'g3qd11<dy>',
        'g3qd12<dx>',
        'g3qd12<dy>',
        'g3qd21<dx>',
        'g3qd21<dy>',
        'g3qd22<dx>',
        'g3qd22<dy>',
    ])

    dx_h3 = parse_errors([
        'h3qd11<dx>',
        'h3qd11<dy>',
        'h3qd12<dx>',
        'h3qd12<dy>',
        'h3qd21<dx>',
        'h3qd21<dy>',
        'h3qd22<dx>',
        'h3qd22<dy>',
    ])

    dx_b3 = parse_errors([
        'b3qd11<dx>',
        'b3qd11<dy>',
        'b3qd12<dx>',
        'b3qd12<dy>',
        'b3qk2<dx>',
        'b3qk2<dy>',
    ])

    rot_ang = parse_errors(['xant_rot_angle', 'yant_rot_angle'])

    mon_first = ['g3dg3g']
    mon_second = ['g3dg5g', 't3dg2g', 't3dg1g', 't3df1']
    mon_gantry = mon_first + mon_second

    # Nelder-Mead
    # CG
    # BFGS
    # L-BFGS-B
    # TNC
    options = dict(
        algorithm='minimize',
        #method='Nelder-Mead',
        mode='xy',
        iterations=30,
        delta=1e-5,
        bounds=Bounds(-0.001, 0.001),
        tol=1e-10)


    ana.fit(init_tw + dx_b3 + dx_h3 + dx_first, ana.monitors, save_to='plots/2-fit.txt', **options)
#   ana.plot_monitors(mon_gantry, save_to='plots/2-fit')
#   ana.plot_steerers(save_to='plots/2-fit')

#   ana.fit(kl_second, mon_second, save_to='plots/3-fit.txt', **options)

#   ana.fit(init_tw + kl_first + kl_second, mon_gantry, **options)
#   ana.fit(init_tw + kl_first + kl_second + dx_first, ana.monitors,
#           **options)

#   ana.info("after-fit-g3dg3g")
#   ana.plot_monitors(save_to='plots/2-fit')
#   ana.plot_steerers(save_to='plots/2-fit')

    ana.plot_monitors(mon_gantry, save_to='plots/3-final')
#   ana.plot_steerers(save_to='plots/3-final')


# monitors:
#   - h1dg1g
#   - h1dg2g
#   - h2dg2g
#   - h3dg3g
#   - b3dg2g
#   - g3dg3g
#   - g3dg5g
#   - t3dg2g
#   - t3dg1g
#   - t3df1

# quads:
#   kL_S0QG1F
#   kL_S0QG4D
#   kl_h1qd11*
#   kl_h1qd12*
#   kl_h2qt11*
#   kl_h2qt12*
#   kl_h2qt13*
#   kl_h3qd11*
#   kl_h3qd12*
#   kl_h3qd21*
#   kl_h3qd22*
#   kl_b3qd11*
#   kl_b3qd12*
#   kl_b3qk2*
#   kl_g3qd11*
#   kl_g3qd12*
#   kl_g3qd21*
#   kl_g3qd22*
#   kl_g3qd31*
#   kl_g3qd32*
#   kl_efg_g3qd41*
#   kl_efg_g3qd42*


# other:
#   g3qd11->tilt: 0.01
#   g3qd12->tilt: 0.01
#   g3qd21->tilt: 0.01
#   g3qd22->tilt: 0.01
#   g3qd31->tilt: 0.01
#   g3qd32->tilt: 0.01

#   scale_error: 0.1
#   gant_rot->angle: 0.01
#   g3qd11->k1: 0.01
#   g3qd12->k1: 0.01
#   dax_G3MU1: 0.0001
#   g3mu1->angle: 0.01
#   g3mu1->e1: 0.01
#   g3mu1->e2: 0.01
#   g3mu1->k0: 0.01
#   g3mu1->hgap: 0.1
#   g3mu1->fint: 0.1
#   g3qd21->k1: 0.01
#   g3ms1v->kick: 0.0001
#   g3qd22->k1: 0.01

#   g3qd11/g3qd22<dx>
#   g3qd11/g3qd22<dy>
#   g3qd11<dx>
#   g3qd11<dy>
#   g3qd12<dx>
#   g3qd12<dy>
#   g3mu1<dx>
#   g3mu1<dy>
#   g3qd21<dx>
#   g3qd21<dy>
#   g3ms1v<dx>
#   g3ms1v<dy>
#   g3qd22<dx>
#   g3qd22<dy>

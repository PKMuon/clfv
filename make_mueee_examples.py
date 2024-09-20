#!/usr/bin/env python3

import numpy as np
from MG5Run import pool
from make_mueee_point import mueee_st_card, zp_width_spline, zp_masses, muon_energies

args = []
for zp_mass, muon_energy in ((0.11, 0.93), (0.15, 11.1), (0.20, 28.2)):
    args.append({
        'workdir': f'mueee_example_{zp_mass}_{muon_energy}',
        'nevent': 100000, 'zp_mass': zp_mass, 'zp_width': zp_width_spline.compute(zp_mass)[0],
        'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
    })
pool.map(mueee_st_card.run, args)

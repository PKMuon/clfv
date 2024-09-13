#!/usr/bin/env python3

import numpy as np
from MG5Run import pool
from make_muemumu_point import muemumu_st_card, zp_width_spline, zp_masses, muon_energies

args = []
for zp_mass, muon_energy in ((0.22, 33.6), (0.25, 50.2), (0.30, 77.2)):
    args.append({
        'workdir': f'muemumu_example_{zp_mass}_{muon_energy}',
        'nevent': 100000, 'zp_mass': zp_mass, 'zp_width': zp_width_spline.compute(zp_mass)[0],
        'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
    })
pool.map(muemumu_st_card.run, args)

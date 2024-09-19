#!/usr/bin/env python3

import numpy as np
from MG5Run import pool
from make_mupep_point import mupep_st_card, zp_width_spline, zp_masses, muon_energies

args = []
for zp_mass, muon_energy in ((0.22, 0.2), (0.22, 2), (0.22, 20)):
    args.append({
        'workdir': f'mupep_example_{zp_mass}_{muon_energy}',
        'nevent': 100000, 'zp_mass': zp_mass, 'zp_width': zp_width_spline.compute(zp_mass)[0],
        'proton_energy': 0.938, 'muon_energy': muon_energy,
    })
pool.map(mupep_st_card.run, args)

#!/usr/bin/env python3

import numpy as np
from make_muemumu_point import make_muemumu_ki_point, zp_masses, muon_energies
from MG5Run import pool
zp_masses = np.linspace(0.01, 1.00, 100)

args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
pool.map(make_muemumu_ki_point, args)

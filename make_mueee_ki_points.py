#!/usr/bin/env python3

import numpy as np
from make_mueee_point import make_mueee_ki_point, zp_masses, muon_energies
from MG5Run import pool

args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
pool.map(make_mueee_ki_point, args)
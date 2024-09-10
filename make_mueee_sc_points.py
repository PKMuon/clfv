#!/usr/bin/env python3

import numpy as np
from make_mueee_point import make_mueee_sc_point, zp_masses, muon_energies
from MG5Run import pool

args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
xses = np.array(pool.map(make_mueee_sc_point, args))
xses = xses.reshape(len(zp_masses), len(muon_energies), -1)
np.savez('data/mueee_sc_points.npz', zp_masses=zp_masses, muon_energies=muon_energies, xses=xses)

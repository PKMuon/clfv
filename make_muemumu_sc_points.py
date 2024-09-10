#!/usr/bin/env python3

import numpy as np
from make_muemumu_point import make_muemumu_sc_point, zp_masses, muon_energies
from MG5Run import pool
muon_energies = np.logspace(1.52, 2, 49)

args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
xses = np.array(pool.map(make_muemumu_sc_point, args))
xses = xses.reshape(len(zp_masses), len(muon_energies), -1)
np.savez('data/muemumu_sc_points.npz', zp_masses=zp_masses, muon_energies=muon_energies, xses=xses)

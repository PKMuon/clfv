#!/usr/bin/env python3

import numpy as np
import multiprocessing
from make_mueee_point import make_mueee_tc_point, zp_masses, muon_energies

pool = multiprocessing.Pool()
args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
xses = np.array(pool.map(make_mueee_tc_point, args))
xses = xses.reshape(len(zp_masses), len(muon_energies), -1)
np.savez('data/mueee_tc_points.npz', zp_masses=zp_masses, muon_energies=muon_energies, xses=xses)

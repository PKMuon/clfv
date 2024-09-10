#!/usr/bin/env python3

import numpy as np
import multiprocessing
from make_muemumu_point import make_muemumu_st_point, zp_masses, muon_energies

pool = multiprocessing.Pool()
args = [ ]
for zp_mass in zp_masses:
    for muon_energy in muon_energies:
        args.append((zp_mass, muon_energy))
xses = np.array(pool.map(make_muemumu_st_point, args))
xses = xses.reshape(len(zp_masses), len(muon_energies), -1)
np.savez('data/muemumu_st_points.npz', zp_masses=zp_masses, muon_energies=muon_energies, xses=xses)

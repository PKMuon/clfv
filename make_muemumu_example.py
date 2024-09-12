#!/usr/bin/env python3

import numpy as np
from make_muemumu_point import muemumu_st_card, zp_width_spline, zp_masses, muon_energies

zp_mass = 0.22
muon_energy = 33.6
zp_width = zp_width_spline.compute(zp_mass)[0]
muemumu_st_card.run({
    'workdir': 'muemumu_example',
    'nevent': 10000, 'zp_mass': zp_mass, 'zp_width': zp_width,
    'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
})

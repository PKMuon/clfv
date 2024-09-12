#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

for channel in 'sc', 'tc':
    muemumu_points = np.load(f'data/muemumu_{channel}_points.npz')
    zp_masses = muemumu_points['zp_masses']
    muon_energies = muemumu_points['muon_energies']
    xses = muemumu_points['xses']
    for zp_mass, xs in zip(zp_masses, xses):
        prefix = r'$m_{Z^\prime} = $' + f'{zp_mass} GeV, {channel[0]}-channel'
        #plt.plot(muon_energies, xs, '.-', label=[prefix, prefix + ', uncertainty'])
        #plt.plot(muon_energies, xs[:,0], '.-', label=prefix)
        #plt.plot(muon_energies, xs[:,0], label=prefix)
        plt.plot(muon_energies, xs[:,0], '-' if channel == 'sc' else '--', label=prefix)
plt.xlabel('Muon energy [GeV]')
plt.ylabel('Cross section [pb]')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig('plot_muemumu_points.pdf')

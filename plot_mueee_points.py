#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

for channel in 'sc', 'tc':
    mueee_points = np.load(f'data/mueee_{channel}_points.npz')
    zp_masses = mueee_points['zp_masses']
    muon_energies = mueee_points['muon_energies']
    xses = mueee_points['xses']
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
plt.legend(loc='center left')
plt.grid()
plt.tight_layout()
plt.savefig('plot_mueee_points.pdf')

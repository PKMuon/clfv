#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

mupep_points = np.load(f'data/mupep_st_points.npz')
zp_masses = mupep_points['zp_masses']
muon_energies = mupep_points['muon_energies']
xses = mupep_points['xses']
for i, (zp_mass, xs) in enumerate(zip(zp_masses, xses)):
    prefix = r'$m_{Z^\prime} = $' + f'{zp_mass:.2f} GeV'
    #plt.plot(muon_energies, xs, '.-', label=[prefix, prefix + ', uncertainty'])
    #plt.plot(muon_energies, xs[:,0], '.-', label=prefix)
    plt.plot(muon_energies, xs[:,0] * 10**i, label=prefix + r' ($\times$' + f'{10**i})')
plt.xlabel('Muon energy [GeV]')
plt.ylabel(r'Cross section [pb]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('plot_mupep_points_st.pdf')

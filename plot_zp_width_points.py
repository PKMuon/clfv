#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

zp_width_points = np.load('data/zp_width_points.npz')
zp_masses = zp_width_points['zp_masses']
zp_widths = zp_width_points['zp_widths']
plt.plot(zp_masses, zp_widths, '.-', label=['width', 'uncertainty'])
plt.xlabel('Mass [GeV]')
plt.ylabel('Width [GeV]')
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('plot_zp_width_points.pdf')

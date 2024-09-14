#!/usr/bin/env python3

import re
import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ki_points = glob.glob('data/muemumu_ki_*.root')
ki_point_pattern = re.compile(r'muemumu_ki_([^_/]*)_([^_/]*)\.root$')

# Load kinematics points.
ki_point_index = { }  # [zp_mass][muon_energy]
for ki_point in ki_points:
    ki_point_match = ki_point_pattern.search(ki_point)
    zp_mass, muon_energy = map(float, ki_point_match.groups())
    if zp_mass not in ki_point_index: ki_point_index[zp_mass] = { }
    ki_point_index[zp_mass][muon_energy] = ki_point
for zp_mass in ki_point_index:
    ki_point_index[zp_mass] = sorted(ki_point_index[zp_mass].items())
ki_points = sorted(ki_point_index.items())

# Draw a plot for each (zp_mass, hist_name).
for zp_mass, zp_mass_points in ki_points:
    lc_index = { }  # [hist_name]
    lc_meta = { }  # [hist_name]
    muon_energies = [ ]
    for muon_energy, ki_point in zp_mass_points:
        muon_energies.append(muon_energy)
        print(ki_point)
        with uproot.open(ki_point) as ki_point:
            for hist_name in ki_point.keys():
                y, x = ki_point[hist_name].to_numpy()
                x = np.mean([x[:-1] + x[1:]], axis=0)
                hist_name = re.sub(r';.*$', '', hist_name)
                if hist_name not in lc_index:
                    lc_index[hist_name] = [ ]
                    lc_meta[hist_name] = {
                        'xtitle': ki_point[hist_name].member("fXaxis").member("fTitle"),
                        'ytitle': ki_point[hist_name].member("fYaxis").member("fTitle"),
                    }
                lc_index[hist_name].append(np.column_stack([x, y]))
    muon_energies = np.array(muon_energies)

    for hist_name, lc in lc_index.items():
        colors = plt.cm.viridis(np.linspace(0, 1, len(muon_energies)))
        lc = LineCollection(lc, colors=colors)
        plt.figure()
        plt.gca().add_collection(lc)
        plt.gca().autoscale()
        plt.xlabel(lc_meta[hist_name]['xtitle'])
        plt.ylabel(lc_meta[hist_name]['ytitle'])
        cbar = plt.colorbar(lc)
        tickposes = np.arange(len(muon_energies))[::10]
        if len(tickposes) * 10 != len(muon_energies): tickposes = np.append(tickposes, len(muon_energies) - 1)
        cbar.set_ticks(np.linspace(0, 1, len(tickposes)))
        cbar.set_ticklabels(['%.2e' % f for f in muon_energies[tickposes]])
        plt.tight_layout()
        plt.savefig(hist_name + '_' + str(zp_mass) + '.pdf')
        plt.close()

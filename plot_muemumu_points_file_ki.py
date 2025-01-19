#!/usr/bin/env python3

import ctypes
import re
import glob
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ki_points = glob.glob('data/muemumu_ki_*.root')
ki_point_pattern = re.compile(r'muemumu_ki_([^_/]*)\.root$')

# Load kinematics points.
ki_point_index = { }  # [zp_mass][muon_energy]
for ki_point in ki_points:
    ki_point_match = ki_point_pattern.search(ki_point)
    if not ki_point_match: continue
    zp_mass = float(ki_point_match.group(1))
    ki_point_index[zp_mass] = ki_point
ki_points = sorted(ki_point_index.items())

# Draw a plot for each (zp_mass, hist_name).
for zp_mass, zp_mass_points in ki_points:
    if zp_mass not in [0.22, 0.3]: continue
    lc_index = { }  # [hist_name]
    lc_meta = { }  # [hist_name]
    print(zp_mass_points)
    with ROOT.TFile.Open(zp_mass_points) as zp_mass_points:
        lc_hists = { }  # [hist_name]
        tree = zp_mass_points.Get('ki_points')
        muon_energy_branch = np.empty(1, dtype='float64')
        tree.SetBranchAddress('muon_energy', muon_energy_branch)
        for branch in tree.GetListOfBranches():
            if branch.GetTitle().endswith('/D'): continue
            hist = ROOT.TH1D()
            lc_hists[branch.GetName()] = hist
            branch.SetAddress(ROOT.AddressOf(hist))
        muon_energies = [ ]
        ientry = 0
        while tree.GetEntry(ientry):
            ientry += 1
            muon_energies.append(muon_energy_branch[0])
            for hist_name, hist in lc_hists.items():
                data = []
                for i in range(1, hist.GetNbinsX() + 1): data.append([hist.GetBinCenter(i), hist.GetBinContent(i)])
                if hist_name not in lc_index:
                    lc_index[hist_name] = [ ]
                    lc_meta[hist_name] = {
                        'xtitle': hist.GetXaxis().GetTitle(),
                        'ytitle': hist.GetYaxis().GetTitle(),
                    }
                lc_index[hist_name].append(np.array(data))
        muon_energies = np.array(muon_energies)
    for hist_name, lc in lc_index.items():
        colors = plt.cm.viridis(np.linspace(0, 1, len(muon_energies)))
        lc = LineCollection(lc, colors=colors)
        plt.figure(figsize=(4, 3))
        plt.gca().add_collection(lc)
        plt.gca().autoscale()
        plt.xlabel(r'$\alpha$ [mrad]')
        plt.ylabel(lc_meta[hist_name]['ytitle'])
        cbar = plt.colorbar(lc)
        cbar.set_label(r'$E_\mu$ [TeV]')
        tickposes = np.arange(len(muon_energies))[::20]
        if len(tickposes) * 20 != len(muon_energies): tickposes = np.append(tickposes, len(muon_energies) - 1)
        cbar.set_ticks(np.linspace(0, 1, len(tickposes)))
        cbar.set_ticklabels(['%.3f' % (f / 1000) for f in muon_energies[tickposes]])
        plt.tight_layout()
        plt.savefig(hist_name + '_' + '%.2f' % zp_mass + '.pdf')
        plt.close()

#!/usr/bin/env python3

import re
import glob
import math
import numpy as np
import ROOT

ki_points = glob.glob('data/muemumu_ki_*_*.root')
ki_point_pattern = re.compile(r'muemumu_ki_([^_/]*)_([^_/]*)\.root$')
muemumu_st_points = np.load('data/muemumu_st_points.npz')

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

for zp_mass, zp_mass_points in ki_points:
    i = muemumu_st_points['zp_masses'].searchsorted(zp_mass)

    outfile = ROOT.TFile.Open(f'data/muemumu_ki_{zp_mass:.18e}.root', 'recreate')
    tree = ROOT.TTree('ki_points', 'kinematic points')
    zp_mass_branch = np.array(zp_mass, dtype=np.float64)
    muon_energy_branch = np.empty(1, dtype=np.float64)
    xs_branch = np.empty(1, dtype=np.float64)
    xs_unc_branch = np.empty(1, dtype=np.float64)
    alpha_pos_branch = ROOT.TH1D()
    alpha_neg_branch = ROOT.TH1D()
    tree.Branch('zp_mass', zp_mass_branch, 'zp_mass/D')
    tree.Branch('muon_energy', muon_energy_branch, 'muon_energy/D')
    tree.Branch('xs', xs_branch, 'xs/D')
    tree.Branch('xs_unc', xs_unc_branch, 'xs_unc/D')
    tree.Branch('alpha_pos', alpha_pos_branch)
    tree.Branch('alpha_neg', alpha_neg_branch)

    for muon_energy, ki_point in zp_mass_points:
        j = muemumu_st_points['muon_energies'].searchsorted(muon_energy)
        muon_energy_branch[0] = muon_energy
        xs_branch[0], xs_unc_branch[0] = muemumu_st_points['xses'][i,j]
        file = ROOT.TFile.Open(ki_point)
        alpha_pos = file.Get('alpha_pos')
        alpha_neg = file.Get('alpha_neg')
        tree.SetBranchAddress('alpha_pos', alpha_pos)
        tree.SetBranchAddress('alpha_neg', alpha_neg)
        tree.Fill()
        file.Close()

    outfile.cd()
    tree.Write()
    outfile.Close()

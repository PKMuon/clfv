#!/usr/bin/env python3

import os
import sys
if len(sys.argv) == 1:
    import glob
    for d in glob.glob('muemumu_example_*/'):
        sys.argv = [sys.argv[0], d + 'Events/run_01/unweighted_events.lhe.gz', d + 'Events/run_01/distributions.root']
        exec(open(sys.argv[0]).read())
    sys.exit()
else:
    assert len(sys.argv) == 3
    INFILE = sys.argv[1]
    OUTFILE = sys.argv[2]

import lhereader
import ROOT
import math

if INFILE.endswith('.gz'):
    import gzip
    with gzip.open(INFILE) as gzfile:
        with open(INFILE[:-3], 'wb') as file:
            file.write(gzfile.read())
INFILE = INFILE[:-3]
outfile = ROOT.TFile(OUTFILE, 'RECREATE')
p_pos = ROOT.TH1D('p_pos', 'momentum of #mu^{+};p_{#mu^{+}} [GeV];Events', 51, 0, 0.2)
p_neg = ROOT.TH1D('p_neg', 'momentum of #mu^{-};p_{#mu^{-}} [GeV];Events', 51, 0, 0.2)
alpha_pos = ROOT.TH1D('alpha_pos', '#alpha of #mu^{+};#alpha_{#mu^{+}};Events', 51, 0, math.pi)
alpha_neg = ROOT.TH1D('alpha_neg', '#alpha of #mu^{-};#alpha_{#mu^{-}};Events', 51, 0, math.pi)
phi_pos = ROOT.TH1D('phi_pos', '#phi of #mu^{+};#phi_{#mu^{+}};Events', 51, -math.pi, math.pi)
phi_neg = ROOT.TH1D('phi_neg', '#phi of #mu^{-};#phi_{#mu^{-}};Events', 51, -math.pi, math.pi)
hists = [p_pos, p_neg, alpha_pos, alpha_neg, phi_pos, phi_neg]

reader = lhereader.LHEReader(INFILE)
for event in reader:
    final_state_particles = [p for p in event.particles if p.status == 1]
    assert len(final_state_particles) == 2
    l_pos = [p for p in final_state_particles if p.pdgid == -13][0]
    l_neg = [p for p in final_state_particles if p.pdgid == +13][0]
    p4_pos = ROOT.TLorentzVector(l_pos.px, l_pos.py, l_pos.pz, l_pos.energy)
    p4_neg = ROOT.TLorentzVector(l_neg.px, l_neg.py, l_neg.pz, l_neg.energy)
    p4 = p4_pos + p4_neg
    boost = p4.BoostVector()
    p4_pos.Boost(-boost)
    p4_neg.Boost(-boost)

    p_pos.Fill(p4_pos.P())
    p_neg.Fill(p4_neg.P())
    alpha_pos.Fill(p4_pos.Theta())
    alpha_neg.Fill(p4_neg.Theta())
    phi_pos.Fill(p4_pos.Phi())
    phi_neg.Fill(p4_neg.Phi())

for hist in hists:
    hist.Write()
    canvas = ROOT.TCanvas()
    hist.Draw()
    canvas.SaveAs(os.path.join(os.path.dirname(OUTFILE), hist.GetName() + '.pdf'))
outfile.Close()

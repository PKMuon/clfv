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
p_pos = ROOT.TH1D('log10p_pos', 'logarithmic momentum of #mu^{+};log_{10}p_{#mu^{+}};Events', 100, 1, 2)
p_neg = ROOT.TH1D('log10p_neg', 'logarithmic momentum of #mu^{-};log_{10}p_{#mu^{-}};Events', 100, 1, 2)
theta_pos = ROOT.TH1D('theta_pos', '#theta of #mu^{+};#theta_{#mu^{+}};Events', 100, 0, 0.01)
theta_neg = ROOT.TH1D('theta_neg', '#theta of #mu^{-};#theta_{#mu^{-}};Events', 100, 0, 0.01)
phi_pos = ROOT.TH1D('phi_pos', '#phi of #mu^{+};#phi_{#mu^{+}};Events', 100, -math.pi, math.pi)
phi_neg = ROOT.TH1D('phi_neg', '#phi of #mu^{-};#phi_{#mu^{-}};Events', 100, -math.pi, math.pi)
hists = [p_pos, p_neg, theta_pos, theta_neg, phi_pos, phi_neg]

reader = lhereader.LHEReader(INFILE)
for event in reader:
    final_state_particles = [p for p in event.particles if p.status == 1]
    assert len(final_state_particles) == 2
    mu_pos = [p for p in final_state_particles if p.pdgid == -13][0]
    mu_neg = [p for p in final_state_particles if p.pdgid == +13][0]
    p4_pos = ROOT.TLorentzVector(mu_pos.px, mu_pos.py, mu_pos.pz, mu_pos.energy)
    p4_neg = ROOT.TLorentzVector(mu_neg.px, mu_neg.py, mu_neg.pz, mu_neg.energy)
    p_pos.Fill(math.log10(p4_pos.P()))
    p_neg.Fill(math.log10(p4_neg.P()))
    theta_pos.Fill(p4_pos.Theta())
    theta_neg.Fill(p4_neg.Theta())
    phi_pos.Fill(p4_pos.Phi())
    phi_neg.Fill(p4_neg.Phi())

for hist in hists:
    hist.Write()
    canvas = ROOT.TCanvas()
    hist.Draw()
    canvas.SaveAs(os.path.join(os.path.dirname(OUTFILE), hist.GetName() + '.pdf'))
outfile.Close()

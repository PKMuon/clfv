import MG5Card
import MG5Width
import MG5Units
import numpy as np

zp_masses = np.array([0.22, 0.25, 0.30])
zp_width_spline = MG5Width.MG5Width('data/zp_width_points.npz', masses_key='zp_masses', widths_key='zp_widths')
muemumu_sc_card = MG5Card.MG5Card('cards/muemumu_sc.dat')
muemumu_tc_card = MG5Card.MG5Card('cards/muemumu_tc.dat')
muemumu_st_card = MG5Card.MG5Card('cards/muemumu_st.dat')
muon_energies = np.logspace(-0.9, 2, 300)

def make_muemumu_point(zp_mass, muon_energy, card):
    print('make_muemumu_point:', zp_mass, muon_energy, sep='\t', flush=True)
    zp_width = zp_width_spline.compute(zp_mass)[0]
    xs, xs_unc, unit, nevent = card.run_xs({
        'nevent': 10000, 'zp_mass': zp_mass, 'zp_width': zp_width,
        'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
    })
    unit = MG5Units.units[unit]
    xs *= unit
    xs_unc *= unit
    print('make_muemumu_point:', zp_mass, muon_energy, xs, xs_unc, sep='\t', flush=True)
    return xs, xs_unc

def make_muemumu_sc_point(args):
    return make_muemumu_point(*args, muemumu_sc_card)

def make_muemumu_tc_point(args):
    return make_muemumu_point(*args, muemumu_tc_card)

def make_muemumu_st_point(args):
    return make_muemumu_point(*args, muemumu_st_card)

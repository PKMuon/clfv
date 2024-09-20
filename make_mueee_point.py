import MG5Card
import MG5Width
import MG5Units
import numpy as np
import os

zp_masses = np.array([0.11, 0.15, 0.20])
zp_width_spline = MG5Width.MG5Width('data/zp_width_points.npz', masses_key='zp_masses', widths_key='zp_widths')
mueee_sc_card = MG5Card.MG5Card('cards/mueee_sc.dat')
mueee_tc_card = MG5Card.MG5Card('cards/mueee_tc.dat')
mueee_st_card = MG5Card.MG5Card('cards/mueee_st.dat')
mueee_ki_card = MG5Card.MG5Card('cards/mueee_ki.dat')
muon_energies = np.logspace(-0.97, 3.00, 498)

def make_mueee_point(zp_mass, muon_energy, card, nevent):
    print('make_mueee_point:', zp_mass, muon_energy, sep='\t', flush=True)
    zp_width = zp_width_spline.compute(zp_mass)[0]
    mueee_kinematics = os.path.abspath(f'data/mueee_ki_{zp_mass:.18e}_{muon_energy:.18e}.root')
    if card is mueee_ki_card and os.path.exists(mueee_kinematics): return
    xs, xs_unc, unit, nevent = card.run_xs({
        'nevent': nevent, 'zp_mass': zp_mass, 'zp_width': zp_width,
        'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
        'make_mueee_distributions': os.path.abspath('make_muell_distributions.py'),
        'mueee_kinematics': mueee_kinematics,
    })
    unit = MG5Units.units[unit]
    xs *= unit
    xs_unc *= unit
    print('make_mueee_point:', zp_mass, muon_energy, xs, xs_unc, sep='\t', flush=True)
    return xs, xs_unc

def make_mueee_sc_point(args):
    return make_mueee_point(*args, mueee_sc_card, 10000)

def make_mueee_tc_point(args):
    return make_mueee_point(*args, mueee_tc_card, 10000)

def make_mueee_st_point(args):
    return make_mueee_point(*args, mueee_st_card, 10000)

def make_mueee_ki_point(args):
    return make_mueee_point(*args, mueee_ki_card, 100000)

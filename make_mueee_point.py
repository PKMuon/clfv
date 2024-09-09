import MG5Card
import MG5Width
import MG5Units
import numpy as np

zp_masses = np.array([0.11, 0.15, 0.20])
zp_width_spline = MG5Width.MG5Width('data/zp_width_points.npz', masses_key='zp_masses', widths_key='zp_widths')
mueee_sc_card = MG5Card.MG5Card('cards/mueee_sc.dat')
mueee_tc_card = MG5Card.MG5Card('cards/mueee_tc.dat')
mueee_st_card = MG5Card.MG5Card('cards/mueee_st.dat')
muon_energies = np.logspace(-0.9, 2, 300)

def make_mueee_point(zp_mass, muon_energy, card):
    print('make_mueee_point:', zp_mass, muon_energy, sep='\t')
    zp_width = zp_width_spline.compute(zp_mass)[0]
    xs, xs_unc, unit, nevent = card.run_xs({
        'nevent': 10000, 'zp_mass': zp_mass, 'zp_width': zp_width,
        'electron_energy': 0.511e-3, 'muon_energy': muon_energy,
    })
    unit = MG5Units.units[unit]
    xs *= unit
    xs_unc *= unit
    print('make_mueee_point:', zp_mass, muon_energy, xs, xs_unc, sep='\t')
    return xs, xs_unc

def make_mueee_sc_point(args):
    return make_mueee_point(*args, mueee_sc_card)

def make_mueee_tc_point(args):
    return make_mueee_point(*args, mueee_tc_card)

def make_mueee_st_point(args):
    return make_mueee_point(*args, mueee_st_card)

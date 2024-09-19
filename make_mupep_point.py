import MG5Card
import MG5Width
import MG5Units
import numpy as np
import os

zp_masses = np.array([0.22, 0.25, 0.30])
zp_width_spline = MG5Width.MG5Width('data/zp_width_points.npz', masses_key='zp_masses', widths_key='zp_widths')
mupep_st_card = MG5Card.MG5Card('cards/mupep_st.dat')
muon_energies = np.logspace(-0.99, 3.00, 400)

def make_mupep_point(zp_mass, muon_energy, card, nevent):
    print('make_mupep_point:', zp_mass, muon_energy, sep='\t', flush=True)
    zp_width = zp_width_spline.compute(zp_mass)[0]
    mupep_kinematics = os.path.abspath(f'data/mupep_ki_{zp_mass:.18e}_{muon_energy:.18e}.root')
    xs, xs_unc, unit, nevent = card.run_xs({
        'nevent': nevent, 'zp_mass': zp_mass, 'zp_width': zp_width,
        'proton_energy': 0.938, 'muon_energy': muon_energy,
    })
    unit = MG5Units.units[unit]
    xs *= unit
    xs_unc *= unit
    print('make_mupep_point:', zp_mass, muon_energy, xs, xs_unc, sep='\t', flush=True)
    return xs, xs_unc

def make_mupep_st_point(args):
    return make_mupep_point(*args, mupep_st_card, 10000)

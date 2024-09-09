import MG5Card
import numpy as np

zp_masses = np.linspace(0.1, 1.0, 46)
#zp_masses = np.logspace(-1, +3, 41)
zp_width_card = MG5Card.MG5Card('cards/zp_width.dat')

units = {'eV': 1e-9}
units['keV'] = units[ 'eV'] * 1e3
units['MeV'] = units['keV'] * 1e3
units['GeV'] = units['MeV'] * 1e3
units['TeV'] = units['GeV'] * 1e3

def make_zp_width_point(zp_mass):
    print('make_zp_width_point:', zp_mass)
    width, width_unc, unit, nevent = zp_width_card.run_width({'nevent': 10000, 'zp_mass': zp_mass})
    unit = units[unit]
    width *= unit
    width_unc *= unit
    print('make_zp_width_point:', zp_mass, width, width_unc, sep='\t')
    return width, width_unc

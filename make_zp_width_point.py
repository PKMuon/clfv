import MG5Card
import numpy as np

zp_masses = np.linspace(0.1, 0.2, 51)
zp_width_card = MG5Card.MG5Card('cards/zp_width.dat')

units = {'MeV': 1}
units['GeV'] = units['MeV'] * 1e3
units['TeV'] = units['GeV'] * 1e3
units['keV'] = units['MeV'] * 1e-3
units['eV' ] = units['keV'] * 1e-3

def make_zp_width_point(zp_mass):
    print('make_zp_width_point:', zp_mass)
    width, width_unc, unit, nevent = zp_width_card.run_width({'nevent': 10000, 'zp_mass': zp_mass})
    unit = units[unit]
    width *= unit
    width_unc *= unit
    return width, width_unc

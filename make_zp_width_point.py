import MG5Card
import MG5Units
import numpy as np

zp_masses = np.linspace(0.1, 1.0, 46)
#zp_masses = np.logspace(-1, +3, 41)
zp_width_card = MG5Card.MG5Card('cards/zp_width.dat')

def make_zp_width_point(zp_mass):
    print('make_zp_width_point:', zp_mass, sep='\t', flush=True)
    width, width_unc, unit, nevent = zp_width_card.run_width({'nevent': 10000, 'zp_mass': zp_mass})
    unit = MG5Units.units[unit]
    width *= unit
    width_unc *= unit
    print('make_zp_width_point:', zp_mass, width, width_unc, sep='\t', flush=True)
    return width, width_unc

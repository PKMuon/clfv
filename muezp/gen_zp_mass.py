#!/usr/bin/env python3

import numpy as np

# 10 GeV
zp_masses = np.linspace(3e-3, 0.0400, 100)
with open('mup_zp_mass_10.0GeV_pb_median.txt', 'w') as file:
    print(*map(lambda x: '%.6f' % x, zp_masses), sep='\n', file=file)

# 4 GeV
zp_masses = np.linspace(3e-3, 0.0170, 100)
with open('mup_zp_mass_4.0GeV_pb_median.txt', 'w') as file:
    print(*map(lambda x: '%.6f' % x, zp_masses), sep='\n', file=file)

# 2 GeV
zp_masses = np.linspace(3e-3, 0.0085, 100)
with open('mup_zp_mass_2.0GeV_pb_median.txt', 'w') as file:
    print(*map(lambda x: '%.6f' % x, zp_masses), sep='\n', file=file)

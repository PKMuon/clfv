#!/usr/bin/env python3

import numpy as np
import glob

for file in glob.glob('mup_energy_*GeV_pb_median.txt'):
    with open(file) as f:
        data = np.array(list(map(float, f.read().strip().split())))
    if data[0] >= 1000.0: data /= 1000.0
    with open(file, 'w') as f:
        print(*map(lambda x: '%.4f' % x, data), file=f, sep='\n')

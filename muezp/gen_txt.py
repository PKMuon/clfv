#!/usr/bin/env python3

import os
import numpy as np
import re
import glob

float_pattern = r'(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?'
rootfile_pattern = re.compile('^MueZp_(' + float_pattern + ')GeV_Zp_(' + float_pattern + ')GeV\\.root$')
mup_energy_pattern = re.compile('mup_energy_(' + float_pattern + ')GeV_pb_median.txt')

for file in os.listdir():
    r = mup_energy_pattern.match(file)
    if not r: continue
    mup_energy = r.group(1)
    muon_energies = set(open(file).read().strip().split())
    cnt = 0
    rootfiles = { }  # [zp_mass]
    for rootfile in os.listdir('mg5root'):
        r = rootfile_pattern.match(rootfile)
        if not r: continue
        muon_energy, zp_mass = r.groups()
        if muon_energy not in muon_energies: continue
        cnt += 1
        rootfiles[zp_mass] = [*rootfiles.get(zp_mass, []), os.path.abspath(os.path.join('mg5root', rootfile))]
    print(f'{cnt} files collected!')
    for zp_mass, files in rootfiles.items():
        with open(f'mup_energy_{mup_energy}GeV_Zp_{zp_mass}GeV.txt', 'w') as file:
            print(*sorted(files), sep='\n', file=file)

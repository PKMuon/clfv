#!/usr/bin/env python3

import numpy as np
from make_zp_width_point import zp_masses, make_zp_width_point
from MG5Run import pool

zp_widths = np.array(pool.map(make_zp_width_point, zp_masses))
np.savez('data/zp_width_points.npz', zp_masses=zp_masses, zp_widths=zp_widths)

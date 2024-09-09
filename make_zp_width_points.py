#!/usr/bin/env python3

import numpy as np
import multiprocessing
from make_zp_width_point import zp_masses, make_zp_width_point

pool = multiprocessing.Pool()
zp_widths = np.array(pool.map(make_zp_width_point, zp_masses))
np.savez('data/zp_width_points.npz', zp_masses=zp_masses, zp_widths=zp_widths)

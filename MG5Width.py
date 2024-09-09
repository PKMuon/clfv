import numpy as np
from scipy.interpolate import make_interp_spline
import warnings

class MG5Width:

    def __init__(self, path, masses_key='masses', widths_key='widths'):
        self.path = path
        data = np.load(path)
        self.masses = data[masses_key]
        self.widths = data[widths_key]
        self.spline = make_interp_spline(self.masses, self.widths)

    def compute(self, mass):
        if not self.masses[0] <= mass <= self.masses[-1]:
            warnings.warn(f'mass {mass} outside interpolation range [{self.masses[0]}, {self.masses[-1]}]')
        return self.spline(mass)

"""DES Y3 likelihood using H0_local from CLASS."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from desy3_likelihood import DESY3 as DESY3_BASE


class DESY3(DESY3_BASE):
    def get_requirements(self):
        req = super().get_requirements()
        req.update({"H0": None, "H0_local": None, "Omega_m": None})
        return req

    def logp(self, **params_values):
        H0_phys = float(self.provider.get_param("H0"))
        H0_loc = float(self.provider.get_param("H0_local"))
        Om_phys = float(self.provider.get_param("Omega_m"))

        R = H0_phys / H0_loc
        Om_loc = Om_phys * (R * R)

        orig_get_param = self.provider.get_param

        def get_param_patched(name):
            if name == "H0":
                return H0_loc
            if name == "Omega_m":
                return Om_loc
            return orig_get_param(name)

        self.provider.get_param = get_param_patched
        try:
            return super().logp(**params_values)
        finally:
            self.provider.get_param = orig_get_param

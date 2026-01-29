from cobaya.likelihood import Likelihood

class SH0ES_H0local(Likelihood):
    """SH0ES constraint using H0_local from holographic CLASS."""
    H0_shoes: float = 73.04
    H0_shoes_err: float = 1.04
    
    def initialize(self):
        self.log.info(f"SH0ES H0_local: {self.H0_shoes} +/- {self.H0_shoes_err}")
    
    def get_requirements(self):
        return {'H0_local': None}
    
    def logp(self, **params_values):
        H0_local = self.provider.get_param('H0_local')
        return -0.5 * ((H0_local - self.H0_shoes) / self.H0_shoes_err) ** 2

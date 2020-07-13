
import numpy as np
import scipy as sp
import scipy.io as sio

class sampling_strategy:
    def __init__(self, config):
        self.tested = []
        self.config = config

    # Internal sample function.  To be overridden in subclasses.
    def _sample(self,mask=None, parameters=None):
        pass

    # Mask is array of True/False (or 1/0).  If True, exclude from sampling
    def sample(self, mask=None, parameters=None):
        individuals = self._sample(mask, parameters)
        self.tested.append(individuals)
        return individuals

    # Fatigue measure.
    # Useful for performing measurements related to unfairness of random testing... "Why
    # do they keep picking on me?"
    #     individuals, count = strategy.fatigue()
    # Individuals is the array of ids of people tested at least once.  
    # Count is an array of times each corresponding individual was tested.
    def fatigue(self):
        return np.unique(self.tested, return_counts=True)


#---------------------------------------------------------------------
# uniform
#
# Assumes config parameter of initializer contains:
#     N - number of samples to draw at a time
#     P - number of individuals in simulation
class uniform(sampling_strategy):
    def __init__(self, config):
        super().__init__(config)
        self.rng = range(0,config['P'])

    def _sample(self,mask=None, parameters=None):
        if mask is None:
            return np.random.choice(self.rng, size=self.config['N'])
        else:
            return np.random.choice(self.rng, size=self.config['N'], p=(1-mask)/config['P'])


#---------------------------------------------------------------------
# weighted
#
# Permits a variety of weighted testing strategies (non-varying weights)
# Assumes config parameters of initializer contains:
#     N - number of samples to draw at a time
#     weights - array of total weights for each node.
class weighted(sampling_strategy):
    def __init__(self, config):
        super().__init__(config)
        self.W = np.array(config['weights'])
        self.W = self.W/np.sum(self.W)
        self.rng = range(0,self.W.shape[0])

    def _sample(self, mask=None, parameters=None):
        tmpW = self.W
        if mask is not None:
            tmpW[mask] = 0
            tmpW = tmpW/np.sum(tmpW)

        return np.random.choice(self.rng, size=self.config['N'], p=tmpW)


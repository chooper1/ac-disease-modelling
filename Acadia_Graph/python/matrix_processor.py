import numpy as np
import scipy as sp
import scipy.stats as st

# Simple class for processing sets of matrices according to a spec.  
# The spec is a dictionary, where the keys are keys of matrices in the "matrices" argument, 
# and the values are dictionaries that describe what to do to them.
# Right now, the operations are:
#     'pre_scale': number  --> multiply the matrix by this weight (default: 1.0).  Intended for helping
#                              prescale matrices before resampling.  e.g. poisson(C * 4) * 1/7 would be:
#                                'C': {'pre_scale': 4.0, 'pre_resample': 'poisson', 'weight': 1.0/7.0}
#     'pre_resample': text --> resampling method.  Currently, 'poisson' is only one implemented
#     'weight': number  --> multiply the matrix by this weight (default: 1.0)

class matrix_processor:
    def __init__(self):
        pass

    def process(self, N, matrices, spec):
        result = np.zeros((N,N))

        for mat, proc in spec.items():
            tmp = matrices[mat]

            if 'pre_scale' in proc:
                tmp = tmp * proc['pre_scale']

            if 'pre_resample' in proc:
                if proc['pre_resample'] == 'poisson':
                    tmp = st.poisson.rvs(tmp)

            if 'weight' in proc:
                tmp = tmp * proc['weight']

            result = result + tmp

        return result


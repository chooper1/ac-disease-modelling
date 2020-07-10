
import numpy as np
import scipy as sp
import scipy.stats as st
import pandas as pd

class testing_facility:
    def __init__(self, fn_rate, fp_rate, turnaround_dist, config=None):
        self.fn_rate = fn_rate
        self.fp_rate = fp_rate
        rng = np.arange(0,turnaround_dist.shape[0])
        self.dist = st.rv_discrete(values=(rng,turnaround_dist))
        self.config = config
        self.result_data = None


    # Submit individuals for tests on given day.  true_results is a true false array
    # of the true test result.
    def submit(self, day, individuals, true_results):
        # Calculate test results
        fp = np.random.uniform(size=len(individuals)) < fp_rate
        fn = np.random.uniform(size=len(individuals)) < fn_rate

        # TODO:  add measure for whether test can work for individual
        # Configuration parameter for test effectiveness by number of days infected?
        # Also need to hand in number of days each individual is infected?
        test_results = np.logical_or(
            np.logical_and(np.logical_not(true_results), fp),
            np.logical_and(true_results, np.logical_not(fn))
            )

        return_times = self.dist.rvs(size=len(individuals)) + day

        tmp = pd.DataFrame({
            'ind':         individuals,
            'submitted':   day,
            'returned':    return_times,
            'true_result': true_results,
            'result':      test_results
            })

        if self.result_data is None:
            self.result_data = tmp
        else:
            self.result_data.append(tmp)


    # Get test results released on given day.
    # Returns two arrays, the set of individuals, and their test results
    def results(self, day):
        applicable = self.result_data.where(self.result_data['returned'] == day)
        return applicable['ind'], applicable['result']


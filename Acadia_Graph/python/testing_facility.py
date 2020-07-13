
import numpy as np
import scipy as sp
import scipy.stats as st
import pandas as pd

class testing_facility:
    def __init__(self, fn_rate, fp_rate, turnaround_dist, config=None):
        self.fn_rate = fn_rate
        self.fp_rate = fp_rate
        rng = np.arange(0,len(turnaround_dist))
        self.dist = st.rv_discrete(values=(rng,turnaround_dist))
        self.config = config
        self.result_data = pd.DataFrame({
            'ind':         [-1],
            'submitted':   [-1],
            'returned':    [-1],
            'true_result': [False],
            'result':      [False]
            })

        types = {
            'ind':         int,
            'submitted':   int,
            'returned':    int,
            'true_result': bool,
            'result':      bool
            }

        self.result_data = self.result_data.astype(types)


    # Submit individuals for tests on given day.  true_results is a true false array
    # of the true test result.
    def submit(self, day, individuals, true_results):
        # Calculate test results
        fp = np.random.uniform(size=len(individuals)) < self.fp_rate
        fn = np.random.uniform(size=len(individuals)) < self.fn_rate

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

        self.result_data = self.result_data.append(tmp, ignore_index=True)


    # Get test results released on given day.
    # Returns two arrays, the set of individuals, and their test results
    def results(self, day):
        applicable = self.result_data[self.result_data['returned'] == day]
        return applicable['ind'].values, applicable['result'].values


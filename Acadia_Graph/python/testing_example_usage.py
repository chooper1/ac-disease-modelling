
import sampling_strategy
import numpy as np
import scipy.stats as st

def main():

    print("---------- Uniform random testing -------------")
    strategy = sampling_strategy.uniform({'N': 30, 'P': 3934})

    print("Doing 10 rounds of testing")
    for i in range(0,10):
        print(strategy.sample())

    print("Fatigue Statistics")
    print(strategy.fatigue())
    print("\n\n")

    print("---------- Weighted random testing -------------")
    # Here, we could do something like:
    #     weights = np.sum(contacts,axis=1)
    weights = st.norm.pdf(range(0,3934),loc=1967,scale=100)
    strategy = sampling_strategy.weighted({'N': 30, 'weights': weights})

    print("Doing 10 rounds of testing")
    for i in range(0,10):
        print(strategy.sample())

    print("Fatigue Statistics")
    print(strategy.fatigue())

if __name__ == "__main__":
    main()


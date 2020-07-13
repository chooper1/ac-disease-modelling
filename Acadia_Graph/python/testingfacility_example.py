
import testing_facility as tf
import numpy as np

def main():
    facility = tf.testing_facility(0.03, 0.02, [0.0, 0.5, 0.45, 0.04, 0.01])

    print(facility.results(0))

    real_statuses = np.array([True, False, False, True, False, False, False, False,
                              True, False, False, False, False, False, True])

    sample = np.array([1, 3, 9, 12], dtype=int)
    facility.submit(0, sample, real_statuses[sample])
    print(facility.results(0))
    print(facility.results(1))
    print(facility.results(2))
    print(facility.results(3))
    print(facility.results(4))

if __name__ == "__main__":
    main()

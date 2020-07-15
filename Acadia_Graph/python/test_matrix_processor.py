
import numpy as np
from matrix_processor import matrix_processor

def main():
    m1 = np.diag([1, 2, 3])
    m2 = np.array([[1,0,0], [1, 1, 1], [2, 0, 0]])
    m3 = np.ones((3,3))

    matrices = {
        'x': m1,
        'y': m2,
        'z': m3
        }

    proc = matrix_processor()
    m4 = proc.process(m1.shape[0], matrices, {
        'x': {'weight': 1.5},
        'y': {'weight': 1},
        'z': {'weight': 2}
        })

    print(m4)

if __name__ == "__main__":
    main()


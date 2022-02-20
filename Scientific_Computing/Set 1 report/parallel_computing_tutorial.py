import multiprocessing as mp
import numpy as np
from time import time

# Prepare data
np.random.RandomState(100)
data = np.random.randint(0, 10, size=[200000, 5])

# Solution Without Paralleization


def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count


results = np.zeros(1000)

for row in data:
    results[row] = howmany_within_range(row, minimum=4, maximum=8)

print(results[:10])

# Parallelizing using Pool.apply()

if __name__ == '__main__':
    # Step 1: Init multiprocessing.Pool()
    pool = mp.Pool(mp.cpu_count())

    # Step 2: `pool.apply` the `howmany_within_range()`
    results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

    # Don't forget to close
    pool.close()

print(results[:, 10])




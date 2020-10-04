import numpy as np
from decimal import Decimal
import time


if __name__ == "__main__":

    t0 = time.time()
    for i in range(1000):
        bitscore = -1400
        float_two = Decimal(2)

        x = np.power(float_two, bitscore) * 1300 * 1300
        alignment_score = -np.log10(x)

    t1 = time.time()
    print(f"1000 Alignment scores: {alignment_score} calculated in {round(t1-t0,4)} seconds")



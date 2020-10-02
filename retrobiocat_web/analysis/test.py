import numpy as np
from decimal import Decimal


if __name__ == "__main__":
    bitscore = Decimal(-1400)
    float_two = Decimal(2)

    x = np.power(float_two, bitscore) * 1300 * 1300
    alignment_score = -np.log10(x)
    print(alignment_score)



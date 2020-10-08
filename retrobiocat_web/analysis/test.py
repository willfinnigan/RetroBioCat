import numpy as np

def average_colour(list_colours):
    return list(np.mean(list_colours, axis=0))


if __name__ == "__main__":
    col_1 = [50, 50, 50]
    col_2 = [25, 50, 75]

    avg = average_colour([col_1, col_2])
    print(avg)
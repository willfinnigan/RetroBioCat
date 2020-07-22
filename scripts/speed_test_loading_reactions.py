from retrobiocat_web.retro.rdchiral.main import rdchiralReaction
from rdkit.Chem import AllChem
import time

def time_function(func, args, num=10):
    times = []
    for i in range(num):
        t0 = time.time()
        func(*args)
        t1 = time.time()
        times.append(t1-t0)

    avg_time = sum(times) / len(times)
    avg_time = round(avg_time, 5)
    print(f"Avg time for function = {avg_time} seconds")
    print()
    return avg_time



if __name__ == '__main__':
    test_smarts = "[#6:1][#6H1:2]=[#8:4]>>[#6:1][#6H2:2][#8H1:4]"

    print("Timing rdkit load reaction from smarts..")
    time_function(AllChem.ReactionFromSmarts, [test_smarts])

    print("Timing rdchiral load reaction from smarts..")
    time_function(rdchiralReaction, [test_smarts])
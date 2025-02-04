import random
import time
# from matplotlib import pyplot as plt


def create_inputs():
    n = [5, 10, 100, 500, 1000, 2000, 4000, 5000]
    characters = ['A', 'C', 'G', 'T']
    sequences = []
    for i in range(n[-1]):
        seq = ''.join(random.choices(characters, k=random.choice(n)))
        sequences.append(seq)
    if len(sequences) % 2 != 0:
        sequences.append(''.join(random.choices(characters, k=random.choice(n))))
    if len(sequences) % 2 != 0:
        sequences.append(''.join(random.choices(characters, k=random.choice(n))))
    with open('gen_inputs.txt', 'w') as f:
        for i in range(0, len(sequences) - 1, 2):
            f.write(sequences[i] + ',' + sequences[i+1] + '\n')
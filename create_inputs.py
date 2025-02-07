import random
import time
from solution import sequence_alignment, read_cost_matrix
import matplotlib.pyplot as plt

#generate random sequences
def generate_random_sequences(lengths, num_pairs=10):
    characters = ['A', 'C', 'G', 'T']
    sequence_pairs = {length: [] for length in lengths}

    for length in lengths:
        for _ in range(num_pairs):
            seq1 = ''.join(random.choices(characters, k=length))
            seq2 = ''.join(random.choices(characters, k=length))
            sequence_pairs[length].append((seq1, seq2))

    return sequence_pairs

def create_inputs():
    lengths = [500, 1000, 2000, 4000, 5000]
    sequence_pairs = generate_random_sequences(lengths)

    with open('gen_inputs.txt', 'w') as file:
        for length in lengths:
            for seq1, seq2 in sequence_pairs[length]:
                file.write(f"{seq1},{seq2}\n")

def measure_runtime(sequence_pairs, loss_matrix):
    runtimes = {length: [] for length in sequence_pairs.keys()}

    for length, pairs in sequence_pairs.items():
        for seq1, seq2 in pairs:
            start_time = time.time()
            sequence_alignment(seq1, seq2, loss_matrix)
            end_time = time.time()
            runtimes[length].append(end_time - start_time)

    average_runtimes = {length: sum(times) / len(times) for length, times in runtimes.items()}
    return average_runtimes

def plot_results(input_sizes, results):
    plt.figure(figsize=(10, 6))
    for name in results:
        plt.plot(input_sizes, results[name], label=name)

    plt.xscale('log')  # Adjust if needed
    plt.yscale('log')  # Adjust if needed
    plt.xlabel('Input Size')
    plt.ylabel('Average Runtime (seconds)')
    plt.title('Algorithm Runtime Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    create_inputs()
    sequence_pairs = generate_random_sequences([500, 1000, 2000, 4000, 5000])
    loss_matrix = read_cost_matrix('imp2cost.txt')
    average_runtimes = measure_runtime(sequence_pairs, loss_matrix)

    input_sizes = list(average_runtimes.keys())
    results = {'Sequence Alignment Algorithm': list(average_runtimes.values())}

    plot_results(input_sizes, results)

    for length, avg_time in average_runtimes.items():
        print(f"Average runtime for sequences of length {length}: {avg_time:.4f} seconds")

if __name__ == "__main__":
    main()

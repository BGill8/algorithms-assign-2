import sys
import numpy as np

# Read cost matrix from file
def read_cost_matrix(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except IOError:
        print(f"Error: Could not read file {filename}")
        sys.exit(1)

    loss_matrix = {}
    symbols = lines[0].strip().split(',')[1:]  # First row: column headers (Y symbols)

    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) != len(symbols) + 1:
            print(f"Error: Invalid format in cost matrix file {filename}")
            sys.exit(1)
        x_symbol = parts[0]  # Row header (X symbol)
        costs = list(map(int, parts[1:]))

        for y_symbol, cost in zip(symbols, costs):
            loss_matrix[(x_symbol, y_symbol)] = cost

    return loss_matrix

# Needleman-Wunsch algorithm for sequence alignment
def needleman_wunsch(seq1, seq2, loss_matrix):
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n+1, m+1), dtype=int)

    # Initialize DP table
    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + loss_matrix[(seq1[i-1], '-')]
    for j in range(1, m+1):
        dp[0][j] = dp[0][j-1] + loss_matrix[('-', seq2[j-1])]

    # Fill DP table
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = dp[i-1][j-1] + loss_matrix[(seq1[i-1], seq2[j-1])]
            delete = dp[i-1][j] + loss_matrix[(seq1[i-1], '-')]
            insert = dp[i][j-1] + loss_matrix[('-', seq2[j-1])]
            dp[i][j] = min(match, delete, insert)

    # Traceback
    i, j = n, m
    aligned_seq1, aligned_seq2 = '', ''
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + loss_matrix[(seq1[i-1], seq2[j-1])]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + loss_matrix[(seq1[i-1], '-')]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, dp[n][m]

# Read multiple sequence pairs from file
def read_sequences(filename):
    sequences = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) == 2:
                    sequences.append((parts[0].strip(), parts[1].strip()))
                else:
                    print(f"Error: Invalid format in sequence file {filename}")
                    sys.exit(1)
    except IOError:
        print(f"Error: Could not read file {filename}")
        sys.exit(1)
    return sequences

def main():

    seq_file = "imp2input.txt"
    cost_file = "imp2cost.txt"
    output_file = "imp2output.txt"

    # Read sequences
    sequence_pairs = read_sequences(seq_file)

    # Read cost matrix
    loss_matrix = read_cost_matrix(cost_file)

    # Process each sequence pair
    try:
        with open(output_file, 'w') as f:
            for seq1, seq2 in sequence_pairs:
                aligned_seq1, aligned_seq2, cost = needleman_wunsch(seq1, seq2, loss_matrix)
                f.write(f"{aligned_seq1},{aligned_seq2}:{cost}\n")
    except IOError:
        print(f"Error: Could not write to file {output_file}")
        sys.exit(1)

    print("Alignment complete. Results saved to", output_file)

if __name__ == "__main__":
    main()

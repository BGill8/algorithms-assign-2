import sys
import numpy as np

# Read cost matrix from file
def read_cost_matrix(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    loss_matrix = {}
    symbols = lines[0].strip().split(',')[1:]  # First row: column headers (Y symbols)

    for line in lines[1:]:
        parts = line.strip().split(',')
        x_symbol = parts[0]  # Row header (X symbol)
        costs = list(map(int, parts[1:]))

        for y_symbol, cost in zip(symbols, costs):
            loss_matrix[(x_symbol, y_symbol)] = cost

    return loss_matrix

# Needleman-Wunsch algorithm for sequence alignment
def needleman_wunsch(seq1, seq2, loss_matrix, gap_penalty=1):
    m, n = len(seq1), len(seq2)

    # Initialize DP table
    dp = np.zeros((m + 1, n + 1))
    traceback = np.zeros((m + 1, n + 1), dtype=int)

    # Fill DP table
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] + gap_penalty
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j-1] + gap_penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + loss_matrix.get((seq1[i-1], seq2[j-1]), float('inf'))
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty

            dp[i][j] = min(match, delete, insert)

            if dp[i][j] == match:
                traceback[i][j] = 1  # Diagonal (match/mismatch)
            elif dp[i][j] == delete:
                traceback[i][j] = 2  # Up (gap in seq2)
            else:
                traceback[i][j] = 3  # Left (gap in seq1)

    # Traceback to find alignment
    aligned_seq1, aligned_seq2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == 1:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, dp[m][n]

# Read multiple sequence pairs from file
def read_sequences(filename):
    sequences = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                sequences.append((parts[0].strip(), parts[1].strip()))
    return sequences

def main():
    if len(sys.argv) != 4:
        print("Usage: python align.py <seq_file> <cost_matrix_file> <output_file>")
        sys.exit(1)

    seq_file, cost_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

    # Read sequences
    sequence_pairs = read_sequences(seq_file)

    # Read cost matrix
    loss_matrix = read_cost_matrix(cost_file)

    # Process each sequence pair
    with open(output_file, 'w') as f:
        for seq1, seq2 in sequence_pairs:
            aligned_seq1, aligned_seq2, cost = needleman_wunsch(seq1, seq2, loss_matrix)
            f.write(f"{aligned_seq1},{aligned_seq2}: {cost}\n")

    print("Alignment complete. Results saved to", output_file)

if __name__ == "__main__":
    main()

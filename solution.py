import numpy as np


def parse_cost_matrix(filename="imp2cost.txt"):
    with open(filename, "r") as f:
        lines = [line.strip().split() for line in f.readlines()]

    alphabet = lines[0][1:]  # First row, excluding the first column
    cost_matrix = {}

    for row in lines[1:]:  # Skip the first row (headers)
        char = row[0]  # Row identifier (A, T, G, C, etc.)
        cost_matrix[char] = {alphabet[i]: int(row[i + 1]) for i in range(len(alphabet))}

    return alphabet, cost_matrix


def align_sequences(seq1, seq2, alphabet, cost_matrix):
    n = len(seq1)
    m = len(seq2)

    # Initialize the matrix with zeros
    matrix = np.zeros((n + 1, m + 1))
    #add weights to matricies
    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + cost_matrix[seq1[i - 1]]["-"]
    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + cost_matrix["-"][seq2[j - 1]]

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            matrix[i][j] = min(
                matrix[i - 1][j] + cost_matrix[seq1[i - 1]]["-"],
                matrix[i][j - 1] + cost_matrix["-"][seq2[j - 1]],
                matrix[i - 1][j - 1] + cost_matrix[seq1[i - 1]][seq2[j - 1]]
            )



def parse_sequences(filename="imp2input.txt"):
    with open(filename, "r") as f:
        sequences = [line.strip().split(",") for line in f.readlines()]
    return sequences



def main():
    alphabet, cost_matrix = parse_cost_matrix()
    sequences = parse_sequences()

    for seq1, seq2 in sequences:
        align_sequences(seq1, seq2, alphabet, cost_matrix)

    


if __name__ == "__main__":
    alphabet, cost_matrix = parse_cost_matrix()
    sequences = parse_sequences()

    print("Alphabet:", alphabet)
    print("Cost Matrix:", cost_matrix)
    print("Sequences:", sequences)

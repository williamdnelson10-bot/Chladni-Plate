import numpy as np


# Same values used for the eigenvalue check in the paper
N = 10
mu = 0.33

h = 2 / N
size = N + 1


# Coefficients from the appendix
b1 = -4 * (3 + mu) * (1 - mu)
b2 = 8 * (1 - mu)
b3 = 2 * (2 - mu)
b4 = -4 * (3 - mu)
b5 = -4 * (2 + mu) * (1 - mu)
b6 = 2 * (1 - mu**2)
b7 = 15 - 8 * mu - 5 * mu**2
b8 = 2 * (8 - 4 * mu - 3 * mu**2)


def make_b1():
    matrix = np.zeros((size, size))

    matrix[0, 0] = b1
    matrix[0, 1] = b2

    for i in range(1, size - 1):
        matrix[i, i - 1] = b3
        matrix[i, i] = b4
        matrix[i, i + 1] = b3

    matrix[-1, -2] = b2
    matrix[-1, -1] = b1

    return matrix


def make_b2():
    matrix = np.zeros((size, size))

    matrix[0, 0] = b1 / 2
    matrix[0, 1] = b3

    matrix[1, 0] = b2 / 4
    matrix[1, 1] = b4 / 2
    matrix[1, 2] = b3 / 2

    for i in range(2, size - 2):
        matrix[i, i - 1] = b3 / 2
        matrix[i, i] = b4 / 2
        matrix[i, i + 1] = b3 / 2

    matrix[-2, -3] = b3 / 2
    matrix[-2, -2] = b4 / 2
    matrix[-2, -1] = b2 / 4

    matrix[-1, -2] = b3
    matrix[-1, -1] = b1 / 2

    return matrix


def make_b3():
    matrix = np.zeros((size, size))

    matrix[0, 0] = b5
    matrix[0, 1] = b3

    matrix[1, 0] = b3 / 2
    matrix[1, 1] = -8
    matrix[1, 2] = 2

    for i in range(2, size - 2):
        matrix[i, i - 1] = 2
        matrix[i, i] = -8
        matrix[i, i + 1] = 2

    matrix[-2, -3] = 2
    matrix[-2, -2] = -8
    matrix[-2, -1] = b3 / 2

    matrix[-1, -2] = b3
    matrix[-1, -1] = b5

    return matrix


def make_d1():
    matrix = np.zeros((size, size))

    matrix[0, 0] = -b1
    matrix[0, 1] = b1
    matrix[0, 2] = b6

    matrix[1, 0] = b1 / 2
    matrix[1, 1] = b7
    matrix[1, 2] = b5
    matrix[1, 3] = b6 / 2

    for i in range(2, size - 2):
        matrix[i, i - 2] = b6 / 2
        matrix[i, i - 1] = b5
        matrix[i, i] = b8
        matrix[i, i + 1] = b5
        matrix[i, i + 2] = b6 / 2

    matrix[-2, -4] = b6 / 2
    matrix[-2, -3] = b5
    matrix[-2, -2] = b7
    matrix[-2, -1] = b1 / 2

    matrix[-1, -3] = b6
    matrix[-1, -2] = b1
    matrix[-1, -1] = -b1

    return matrix


def make_d2():
    matrix = np.zeros((size, size))

    matrix[0, 0] = b7
    matrix[0, 1] = b4
    matrix[0, 2] = 2

    matrix[1, 0] = b4 / 2
    matrix[1, 1] = 18
    matrix[1, 2] = -8
    matrix[1, 3] = 1

    for i in range(2, size - 2):
        matrix[i, i - 2] = 1
        matrix[i, i - 1] = -8
        matrix[i, i] = 19
        matrix[i, i + 1] = -8
        matrix[i, i + 2] = 1

    matrix[-2, -4] = 1
    matrix[-2, -3] = -8
    matrix[-2, -2] = 18
    matrix[-2, -1] = b4 / 2

    matrix[-1, -3] = 2
    matrix[-1, -2] = b4
    matrix[-1, -1] = b7

    return matrix


def make_d3():
    matrix = np.zeros((size, size))

    matrix[0, 0] = b8
    matrix[0, 1] = b4
    matrix[0, 2] = 2

    matrix[1, 0] = b4 / 2
    matrix[1, 1] = 19
    matrix[1, 2] = -8
    matrix[1, 3] = 1

    for i in range(2, size - 2):
        matrix[i, i - 2] = 1
        matrix[i, i - 1] = -8
        matrix[i, i] = 20
        matrix[i, i + 1] = -8
        matrix[i, i + 2] = 1

    matrix[-2, -4] = 1
    matrix[-2, -3] = -8
    matrix[-2, -2] = 19
    matrix[-2, -1] = b4 / 2

    matrix[-1, -3] = 2
    matrix[-1, -2] = b4
    matrix[-1, -1] = b8

    return matrix


# Smaller diagonal matrices from the appendix
I1 = np.zeros((size, size))

I1[0, 0] = b6
I1[-1, -1] = b6

for i in range(1, size - 1):
    I1[i, i] = 2

I2 = I1 / 2


# Build the remaining block matrices
B1 = make_b1()
B2 = make_b2()
B3 = make_b3()

D1 = make_d1()
D2 = make_d2()
D3 = make_d3()


# Full matrix M has (N + 1)^2 rows and columns
matrix_size = size**2
M = np.zeros((matrix_size, matrix_size))


def put_block(block_row, block_column, block):
    row_start = block_row * size
    column_start = block_column * size

    M[
        row_start:row_start + size,
        column_start:column_start + size
    ] = block


# First block row
put_block(0, 0, D1)
put_block(0, 1, B1)
put_block(0, 2, I1)


# Second block row
put_block(1, 0, B2)
put_block(1, 1, D2)
put_block(1, 2, B3)
put_block(1, 3, I2)


# Interior block rows
for row in range(2, size - 2):
    put_block(row, row - 2, I2)
    put_block(row, row - 1, B3)
    put_block(row, row, D3)
    put_block(row, row + 1, B3)
    put_block(row, row + 2, I2)


# Second-last block row
put_block(size - 2, size - 4, I2)
put_block(size - 2, size - 3, B3)
put_block(size - 2, size - 2, D2)
put_block(size - 2, size - 1, B2)


# Last block row
put_block(size - 1, size - 3, I1)
put_block(size - 1, size - 2, B1)
put_block(size - 1, size - 1, D1)


# Calculate the eigenvalues
eigenvalues = np.linalg.eigvals(M)

# Tiny negative values can appear from numerical rounding
positive_eigenvalues = []

for value in eigenvalues:
    if value.real > 1e-9:
        positive_eigenvalues.append(value.real)


# lambda = h^4 * omega^2
frequencies = []

for value in positive_eigenvalues:
    omega = np.sqrt(value) / h**2
    frequencies.append(omega)

frequencies.sort()


# Frequencies D'Alessio reports for N = 10
paper_frequencies = [
    14.0,
    26.5,
    39.7,
    62.6,
    72.2,
]


print("Free-edge eigenvalue test")
print("-------------------------")
print(f"N = {N}")
print(f"h = {h}")
print()


# Find the computed frequency closest to each paper value
for paper_frequency in paper_frequencies:
    closest = min(
        frequencies,
        key=lambda frequency: abs(
            frequency - paper_frequency
        )
    )

    print(
        f"Paper: {paper_frequency:5.1f}   "
        f"Computed: {closest:7.3f}"
    )

import numpy as np
from scipy.sparse import bmat, csr_matrix


# Build the free-edge plate matrix
# N = number of grid intervals along each side
# mu = Poisson's ratio

def build_free_edge_matrix(N, mu):

    size = N + 1

    # Coefficients from D'Alessio's free-edge discretization
    b1 = -4 * (3 + mu) * (1 - mu)
    b2 = 8 * (1 - mu)
    b3 = 2 * (2 - mu)
    b4 = -4 * (3 - mu)
    b5 = -4 * (2 + mu) * (1 - mu)
    b6 = 2 * (1 - mu**2)
    b7 = 15 - 8 * mu - 5 * mu**2
    b8 = 2 * (8 - 4 * mu - 3 * mu**2)


    # Build the smaller matrices used in the full plate matrix

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

        return csr_matrix(matrix)


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

        return csr_matrix(matrix)


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

        return csr_matrix(matrix)


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

        return csr_matrix(matrix)


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

        return csr_matrix(matrix)


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

        return csr_matrix(matrix)


    # Matrices for points two grid spaces apart
    I1 = np.zeros((size, size))

    I1[0, 0] = b6
    I1[-1, -1] = b6

    for i in range(1, size - 1):
        I1[i, i] = 2

    I1 = csr_matrix(I1)
    I2 = I1 / 2


    # Calculate each of the smaller matrices
    B1 = make_b1()
    B2 = make_b2()
    B3 = make_b3()

    D1 = make_d1()
    D2 = make_d2()
    D3 = make_d3()


    # Assemble the full plate matrix from the smaller blocks
    blocks = []

    for row in range(size):
        blocks.append([None] * size)


    # First block row
    blocks[0][0] = D1
    blocks[0][1] = B1
    blocks[0][2] = I1


    # Second block row
    blocks[1][0] = B2
    blocks[1][1] = D2
    blocks[1][2] = B3
    blocks[1][3] = I2


    # Interior block rows
    for row in range(2, size - 2):
        blocks[row][row - 2] = I2
        blocks[row][row - 1] = B3
        blocks[row][row] = D3
        blocks[row][row + 1] = B3
        blocks[row][row + 2] = I2


    # Second-to-last block row
    blocks[size - 2][size - 4] = I2
    blocks[size - 2][size - 3] = B3
    blocks[size - 2][size - 2] = D2
    blocks[size - 2][size - 1] = B2


    # Last block row
    blocks[size - 1][size - 3] = I1
    blocks[size - 1][size - 2] = B1
    blocks[size - 1][size - 1] = D1


    # Combine all blocks into one sparse matrix
    return bmat(blocks, format="csr")

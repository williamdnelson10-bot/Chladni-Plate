import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import bmat, csr_matrix, eye
from scipy.sparse.linalg import spsolve


# -------------------------
# Physical plate parameters
# -------------------------

# Our square aluminum plate
plate_side_length = 17 * 0.0254   # 17 inches, converted to meters
thickness = 0.125 * 0.0254        # 1/8 inch, converted to meters

# Aluminum material properties
youngs_modulus = 69e9             # Pa
density = 2700                    # kg/m^3
mu = 0.33                         # Poisson ratio


# Half side length of the plate
a = plate_side_length / 2


# Plate parameter D
D = np.sqrt(
    youngs_modulus * thickness**2
    / (12 * density * (1 - mu**2))
)


# Dimensionless natural angular frequencies
# found using free_edge_test.py
dimensionless_omega_values = [
    13.966,
    26.455,
    39.658,
    62.570,
    72.208,
]


# Convert the dimensionless natural frequencies
# into physical frequencies for our plate, in Hz
frequencies_to_run = [
    round(
        D * omega / (2 * np.pi * a**2),
        1,
    )
    for omega in dimensionless_omega_values
]


# Match each physical frequency to the
# dimensionless frequency used by the simulation
omega_values = dict(
    zip(
        frequencies_to_run,
        dimensionless_omega_values,
    )
)


# -------------------------
# Numerical grid
# -------------------------

# Number of intervals along each side
# N = 40 to start, D'Alessio used 200
N = 200

# Grid spacing
h = 2 / N

# Number of grid points along one side
size = N + 1

# Total number of grid points
number_of_points = size**2

# Coordinates
x = np.linspace(-1, 1, size)
y = np.linspace(-1, 1, size)

# Center of the plate
center = N // 2

# Position of the center in the flattened vector
center_index = center + center * size


# Make sure the plots folder exists
os.makedirs("plots", exist_ok=True)


def build_free_edge_matrix():
    # Constants from the appendix of the paper
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


    # Diagonal matrices from the appendix
    I1 = np.zeros((size, size))

    I1[0, 0] = b6
    I1[-1, -1] = b6

    for i in range(1, size - 1):
        I1[i, i] = 2

    I1 = csr_matrix(I1)
    I2 = I1 / 2


    # Build the block matrices
    B1 = make_b1()
    B2 = make_b2()
    B3 = make_b3()

    D1 = make_d1()
    D2 = make_d2()
    D3 = make_d3()


    # Empty block layout for the full matrix
    blocks = []

    for row in range(size):
        block_row = [None] * size
        blocks.append(block_row)


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


    # Second-last block row
    blocks[size - 2][size - 4] = I2
    blocks[size - 2][size - 3] = B3
    blocks[size - 2][size - 2] = D2
    blocks[size - 2][size - 1] = B2


    # Last block row
    blocks[size - 1][size - 3] = I1
    blocks[size - 1][size - 2] = B1
    blocks[size - 1][size - 1] = D1


    # Join everything into one sparse matrix
    return bmat(
        blocks,
        format="csr",
    )


print("Building free-edge matrix...")

M = build_free_edge_matrix()

print("Matrix finished.")


# Convert M into the dimensionless spatial plate operator
plate_operator = M / h**4

# Identity matrix used in the harmonic equation
identity = eye(
    number_of_points,
    format="csr",
)


# Point force at the center of the plate
force = np.zeros(number_of_points)

# A point source integrated over the center cell
force[center_index] = 1 / h**2


# Grid used for contour plotting
X, Y = np.meshgrid(x, y)


# Store results so they can also be placed in one combined figure
results = {}


for frequency_hz in frequencies_to_run:
    omega0 = omega_values[frequency_hz]

    print(
        f"Solving {frequency_hz} Hz "
        f"(omega = {omega0})..."
    )


    # For a harmonic response u = F cos(omega t),
    # the spatial amplitude satisfies
    #
    # (L - omega^2 I) F = force
    #
    # where L is the discretized plate operator
    system = (
        plate_operator
        - omega0**2 * identity
    )


    # Solve for the displacement amplitude
    response = spsolve(
        system,
        force,
    )


    # Convert the flattened solution back into the square plate
    u = response.reshape(
        (size, size),
        order="F",
    )

    results[frequency_hz] = u


    # -------------------------
    # Individual nodal plot
    # -------------------------

    fig, ax = plt.subplots(
        figsize=(6, 6)
    )

    ax.contour(
        X,
        Y,
        u.T,
        levels=[0],
        colors="black",
    )

    ax.set_title(
        f"Nodal Lines at {frequency_hz} Hz"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_aspect("equal")

    plt.tight_layout()


    # Save the plot
    fig.savefig(
        f"plots/nodal_{frequency_hz}hz.png",
        dpi=300,
    )

    plt.close(fig)


print("All frequencies finished.")


# -------------------------
# Combined plot
# -------------------------

fig, axes = plt.subplots(
    1,
    len(frequencies_to_run),
    figsize=(18, 4),
)


for ax, frequency_hz in zip(
    axes,
    frequencies_to_run,
):
    u = results[frequency_hz]

    ax.contour(
        X,
        Y,
        u.T,
        levels=[0],
        colors="black",
    )

    ax.set_title(
        f"{frequency_hz} Hz"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_aspect("equal")


fig.suptitle(
    "Free-Edge Plate Nodal Patterns"
)

plt.tight_layout()


# Save all five together
fig.savefig(
    "plots/nodal_all_frequencies.png",
    dpi=300,
)


# Show the combined figure
plt.show()
"""Center-forced free-edge Chladni plate: all recorded experimental frequencies.

Based on simulation/actual_case.py in the Chladni-Plate project.
Run from the project root so plots go under plots/experimental_frequencies/.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import bmat, csr_matrix, eye
from scipy.sparse.linalg import spsolve


# -------------------------
# Physical plate parameters
# -------------------------

plate_side_length = 17 * 0.0254       # 17-inch square plate, meters
thickness = (1 / 16) * 0.0254         # 1/16 inch, meters
youngs_modulus = 69e9                 # Aluminum, Pa
density = 2700                        # Aluminum, kg/m^3
mu = 0.33                             # Poisson's ratio

a = plate_side_length / 2            # Half the side length

# Same plate parameter and frequency scaling as actual_case.py
D = np.sqrt(
    youngs_modulus * thickness**2
    / (12 * density * (1 - mu**2))
)


# -------------------------
# Frequencies to simulate
# -------------------------

# Strings preserve the exact labels used for the videos and images.
# Each group follows the order in the experimental records.
frequency_groups = {
    "1_microphone": [
        "102.93", "190.90", "287.80", "447.18",
        "579.89", "637.15", "684.41", "724.19",
    ],
    "2_simulation_predictions": [
        "127.6", "317.0", "414.2", "506.0", "621.8", "757.2",
    ],
    "3_previous_experiment": [
        "100", "195", "250", "273", "360", "558", "604", "675",
    ],
    "4_harmonic_tests": [
        "75", "150", "225", "300", "375", "450", "525", "600",
    ],
    "5_miscellaneous": [
        "240", "224", "344",
    ],
}

# Run all five groups. To test only one group, change this list.
groups_to_run = list(frequency_groups)


# -------------------------
# Numerical grid
# -------------------------

N = 200                            # Same as our earlier code; try 40 for a fast test
h = 2 / N
size = N + 1
number_of_points = size**2

x = np.linspace(-1, 1, size)
y = np.linspace(-1, 1, size)
X, Y = np.meshgrid(x, y)

center = N // 2
center_index = center + center * size

output_folder = "plots/experimental_frequencies"
os.makedirs(output_folder, exist_ok=True)


# -------------------------
# Build the free-edge matrix
# (same discretization as actual_case.py)
# -------------------------

def build_free_edge_matrix():
    # Constants from the free-edge stencil in the paper's appendix
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

    I1 = np.zeros((size, size))
    I1[0, 0] = b6
    I1[-1, -1] = b6
    for i in range(1, size - 1):
        I1[i, i] = 2
    I1 = csr_matrix(I1)
    I2 = I1 / 2

    B1 = make_b1()
    B2 = make_b2()
    B3 = make_b3()
    D1 = make_d1()
    D2 = make_d2()
    D3 = make_d3()

    blocks = []
    for row in range(size):
        blocks.append([None] * size)

    blocks[0][0] = D1
    blocks[0][1] = B1
    blocks[0][2] = I1

    blocks[1][0] = B2
    blocks[1][1] = D2
    blocks[1][2] = B3
    blocks[1][3] = I2

    for row in range(2, size - 2):
        blocks[row][row - 2] = I2
        blocks[row][row - 1] = B3
        blocks[row][row] = D3
        blocks[row][row + 1] = B3
        blocks[row][row + 2] = I2

    blocks[size - 2][size - 4] = I2
    blocks[size - 2][size - 3] = B3
    blocks[size - 2][size - 2] = D2
    blocks[size - 2][size - 1] = B2

    blocks[size - 1][size - 3] = I1
    blocks[size - 1][size - 2] = B1
    blocks[size - 1][size - 1] = D1

    return bmat(blocks, format="csr")


# -------------------------
# Solve the center-forced response
# -------------------------

print("Building free-edge matrix...")
M = build_free_edge_matrix()
plate_operator = M / h**4
identity = eye(number_of_points, format="csr")

# Approximate a point force at the center grid cell.
force = np.zeros(number_of_points)
force[center_index] = 1 / h**2
print("Matrix finished.")

for group_name in groups_to_run:
    labels = frequency_groups[group_name]
    folder = os.path.join(output_folder, group_name)
    os.makedirs(folder, exist_ok=True)

    results = {}
    print(f"\nRunning {group_name} ({len(labels)} frequencies)")

    for label in labels:
        frequency_hz = float(label)
        omega0 = 2 * np.pi * a**2 * frequency_hz / D

        # Harmonic steady-state equation at the EXACT measured frequency:
        # (L - omega^2 I) F = force
        system = plate_operator - omega0**2 * identity
        response = spsolve(system, force)

        # Convert the vector back to a square plate displacement field.
        u = response.reshape((size, size), order="F")
        results[label] = u
        print(f"  Finished {label} Hz")

        # Individual nodal-line plot (black lines, same as older plots).
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.contour(X, Y, u.T, levels=[0], colors="black")
        ax.set_title(f"Nodal Lines at {label} Hz")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect("equal")
        plt.tight_layout()
        fig.savefig(os.path.join(folder, f"{label} Hz.png"), dpi=300)
        plt.close(fig)

    # One 3-column figure containing every frequency in this group.
    columns = 3
    rows = (len(labels) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, figsize=(11, 3.5 * rows))
    axes = np.array(axes).reshape(-1)

    for index, label in enumerate(labels):
        ax = axes[index]
        ax.contour(X, Y, results[label].T, levels=[0], colors="black")
        ax.set_title(f"{label} Hz")
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect("equal")
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])

    for index in range(len(labels), len(axes)):
        axes[index].axis("off")

    fig.suptitle(group_name.replace("_", " ").title(), fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(folder, "combined.png"), dpi=300)
    plt.close(fig)
    print(f"Saved plots in {folder}")

print("\nAll requested experimental frequencies finished.")

import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import eigs


# -------------------------
# Physical plate parameters
# -------------------------

# Our square aluminum plate
plate_side_length = 17 * 0.0254   # 17 inches, converted to meters
thickness = 0.0625 * 0.0254        # 1/16 inch, converted to meters

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


# -------------------------
# Frequency range
# -------------------------

minimum_frequency = 50
maximum_frequency = 700


# Calculate enough of the lowest eigenvalues
# to get past the 700 Hz range.
number_of_eigenvalues = 40


# -------------------------
# Numerical grid
# -------------------------

# Same grid resolution used in actual_case.py
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


# Folder for the natural mode plots
os.makedirs(
    "plots/natural_modes",
    exist_ok=True,
)


# -------------------------
# Build the free-edge matrix
# -------------------------

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
print("Calculating natural modes...")


# -------------------------
# Calculate eigenvalues and eigenvectors
# -------------------------

# Search near zero to obtain the lowest-frequency modes.
# A tiny nonzero sigma avoids sitting directly on
# the zero-frequency rigid-body modes.
eigenvalues, eigenvectors = eigs(
    M,
    k=number_of_eigenvalues,
    sigma=1e-8,
    which="LM",
)


print("Eigenvalue calculation finished.")


# -------------------------
# Convert to physical frequencies
# -------------------------

modes = []
all_positive_frequencies = []


for index in range(len(eigenvalues)):
    eigenvalue = eigenvalues[index]


    # The eigenvalues we use should be real.
    # Ignore tiny numerical imaginary parts.
    if abs(eigenvalue.imag) > 1e-7:
        continue

    eigenvalue = eigenvalue.real


    # Ignore zero-frequency rigid-body modes
    # and tiny negative numerical values.
    if eigenvalue <= 1e-12:
        continue


    # Convert the eigenvalue into the
    # dimensionless angular natural frequency.
    #
    # lambda = h^4 * omega^2
    omega = np.sqrt(
        eigenvalue
    ) / h**2


    # Convert the dimensionless angular frequency
    # into physical frequency in Hz.
    frequency_hz = (
        D * omega
        / (2 * np.pi * a**2)
    )


    all_positive_frequencies.append(
        frequency_hz
    )


    # Only keep modes inside the range
    # that we tested experimentally.
    if (
        frequency_hz < minimum_frequency
        or frequency_hz > maximum_frequency
    ):
        continue


    # Get the eigenvector for this mode.
    mode_vector = eigenvectors[:, index]


    # Numerical eigenvectors can contain a small
    # complex phase even when the physical mode is real.
    #
    # Rotate the vector so that it can be treated
    # as a real vibration shape.
    largest_index = np.argmax(
        np.abs(mode_vector)
    )

    phase = np.angle(
        mode_vector[largest_index]
    )

    mode_vector = (
        mode_vector
        * np.exp(-1j * phase)
    )

    mode_vector = mode_vector.real


    # Normalize the mode so that its largest
    # displacement anywhere on the plate is 1.
    maximum_amplitude = np.max(
        np.abs(mode_vector)
    )

    mode_vector = (
        mode_vector
        / maximum_amplitude
    )


    # Amount of motion at the exact center
    # relative to the maximum motion in this mode.
    center_response = abs(
        mode_vector[center_index]
    )


    # Convert the flattened vector back into
    # the square plate.
    mode_shape = mode_vector.reshape(
        (size, size),
        order="F",
    )


    modes.append(
        {
            "frequency": frequency_hz,
            "omega": omega,
            "center_response": center_response,
            "shape": mode_shape,
        }
    )


# Sort all modes from lowest to highest frequency.
modes.sort(
    key=lambda mode: mode["frequency"]
)


# -------------------------
# Identify repeated frequencies
# -------------------------

# Modes within this tolerance are treated as
# having essentially the same natural frequency.
#
# IMPORTANT:
# We only GROUP them for labeling.
# We do NOT remove any of their mode shapes.
frequency_tolerance = 0.1


frequency_groups = []


for mode in modes:
    if len(frequency_groups) == 0:
        frequency_groups.append(
            [mode]
        )
        continue


    group_frequency = (
        frequency_groups[-1][0]["frequency"]
    )


    if (
        abs(
            mode["frequency"]
            - group_frequency
        )
        <= frequency_tolerance
    ):
        frequency_groups[-1].append(
            mode
        )

    else:
        frequency_groups.append(
            [mode]
        )


# Add information about repeated frequencies
# to every individual mode.
for group_number, group in enumerate(
    frequency_groups,
    start=1,
):

    for position, mode in enumerate(
        group,
        start=1,
    ):
        mode["group_number"] = group_number
        mode["group_size"] = len(group)
        mode["group_position"] = position


# -------------------------
# Print the results
# -------------------------

print()
print(
    f"Natural modes from "
    f"{minimum_frequency} to "
    f"{maximum_frequency} Hz"
)

print(
    "---------------------------------------------"
)


for mode_number, mode in enumerate(
    modes,
    start=1,
):

    frequency = mode["frequency"]

    center_percent = (
        100 * mode["center_response"]
    )

    group_size = mode["group_size"]
    group_position = mode["group_position"]


    if group_size > 1:

        repeated_text = (
            f"Mode {group_position} of "
            f"{group_size} at this frequency"
        )

    else:

        repeated_text = (
            "Single mode at this frequency"
        )


    print(
        f"{mode_number:2d}: "
        f"{frequency:7.1f} Hz   "
        f"Center = "
        f"{center_percent:5.1f}%   "
        f"{repeated_text}"
    )


# Make sure the eigenvalue calculation
# actually extended above 700 Hz.
if len(all_positive_frequencies) > 0:

    highest_frequency = max(
        all_positive_frequencies
    )

    print()

    print(
        f"Highest frequency calculated: "
        f"{highest_frequency:.1f} Hz"
    )


    if highest_frequency < maximum_frequency:

        print(
            "Increase number_of_eigenvalues "
            "because the calculation did not "
            "reach 700 Hz."
        )


# -------------------------
# Plot every natural mode
# -------------------------

X, Y = np.meshgrid(
    x,
    y,
)


for mode_number, mode in enumerate(
    modes,
    start=1,
):

    frequency = mode["frequency"]
    mode_shape = mode["shape"]

    center_percent = (
        100 * mode["center_response"]
    )

    group_size = mode["group_size"]
    group_position = mode["group_position"]


    fig, ax = plt.subplots(
        figsize=(6, 6)
    )


    ax.contour(
        X,
        Y,
        mode_shape.T,
        levels=[0],
        colors="black",
    )


    # If multiple modes have the same frequency,
    # make that clear in the title.
    if group_size > 1:

        ax.set_title(
            f"{frequency:.1f} Hz "
            f"- Mode "
            f"{group_position} of "
            f"{group_size}\n"
            f"Center = "
            f"{center_percent:.1f}%"
        )

    else:

        ax.set_title(
            f"{frequency:.1f} Hz\n"
            f"Center = "
            f"{center_percent:.1f}%"
        )


    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_aspect("equal")

    plt.tight_layout()


    fig.savefig(
        "plots/natural_modes/"
        f"mode_{mode_number:02d}_"
        f"{frequency:.1f}hz.png",
        dpi=300,
    )


    plt.close(fig)


# -------------------------
# Combined plot
# -------------------------

number_of_modes = len(modes)

columns = 4

rows = (
    number_of_modes
    + columns
    - 1
) // columns


fig, axes = plt.subplots(
    rows,
    columns,
    figsize=(
        columns * 4,
        rows * 4,
    ),
)


# Turn all subplot locations
# into one simple list.
axes = np.array(
    axes
).reshape(-1)


for index in range(number_of_modes):

    ax = axes[index]

    mode = modes[index]

    frequency = mode["frequency"]
    mode_shape = mode["shape"]

    center_percent = (
        100 * mode["center_response"]
    )

    group_size = mode["group_size"]
    group_position = mode["group_position"]


    ax.contour(
        X,
        Y,
        mode_shape.T,
        levels=[0],
        colors="black",
    )


    if group_size > 1:

        ax.set_title(
            f"{frequency:.1f} Hz "
            f"({group_position}/{group_size})\n"
            f"Center = "
            f"{center_percent:.1f}%"
        )

    else:

        ax.set_title(
            f"{frequency:.1f} Hz\n"
            f"Center = "
            f"{center_percent:.1f}%"
        )


    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_aspect("equal")


# Hide any empty subplot spaces.
for index in range(
    number_of_modes,
    len(axes),
):

    axes[index].axis("off")


fig.suptitle(
    "Predicted Natural Modes: "
    "50-700 Hz"
)


plt.tight_layout()


fig.savefig(
    "plots/natural_modes/"
    "natural_modes_50_700hz.png",
    dpi=300,
)


plt.show()
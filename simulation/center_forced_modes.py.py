import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import bmat, csr_matrix, identity
from scipy.sparse.linalg import eigs, spsolve


# -------------------------
# Physical plate parameters
# -------------------------

# Our square aluminum plate
plate_side_length = 17 * 0.0254   # 17 inches in meters
thickness = 0.0625 * 0.0254        # 1/16 inch, converted to meters

# Aluminum material properties
youngs_modulus = 69e9             # Pa
density = 2700                    # kg/m^3
mu = 0.33                         # Poisson ratio


# Half side length
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
maximum_frequency = 1000


# Same idea as natural_modes.py.
# This should calculate enough eigenvalues
# to get beyond 2000 Hz.
number_of_eigenvalues = 100


# Frequencies this close together are
# treated as the same resonance.
frequency_tolerance = 0.1


# The forced equation becomes singular
# exactly at an ideal undamped resonance.
#
# So we solve just slightly above it.
resonance_offset = 0.1


# -------------------------
# Numerical grid
# -------------------------

N = 200

h = 2 / N

size = N + 1

number_of_points = size**2


x = np.linspace(
    -1,
    1,
    size,
)

y = np.linspace(
    -1,
    1,
    size,
)


X, Y = np.meshgrid(
    x,
    y,
)


# Center of the plate
center = N // 2

center_index = (
    center
    + center * size
)


# Folder for plots
os.makedirs(
    "plots/center_forced_modes",
    exist_ok=True,
)


# -------------------------
# Build the free-edge matrix
# -------------------------

def build_free_edge_matrix():

    # Constants from the appendix
    b1 = -4 * (3 + mu) * (1 - mu)
    b2 = 8 * (1 - mu)
    b3 = 2 * (2 - mu)
    b4 = -4 * (3 - mu)
    b5 = -4 * (2 + mu) * (1 - mu)
    b6 = 2 * (1 - mu**2)
    b7 = 15 - 8 * mu - 5 * mu**2
    b8 = 2 * (8 - 4 * mu - 3 * mu**2)


    def make_b1():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = b1
        matrix[0, 1] = b2

        for i in range(
            1,
            size - 1,
        ):

            matrix[i, i - 1] = b3
            matrix[i, i] = b4
            matrix[i, i + 1] = b3

        matrix[-1, -2] = b2
        matrix[-1, -1] = b1

        return csr_matrix(
            matrix
        )


    def make_b2():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = b1 / 2
        matrix[0, 1] = b3

        matrix[1, 0] = b2 / 4
        matrix[1, 1] = b4 / 2
        matrix[1, 2] = b3 / 2

        for i in range(
            2,
            size - 2,
        ):

            matrix[i, i - 1] = b3 / 2
            matrix[i, i] = b4 / 2
            matrix[i, i + 1] = b3 / 2

        matrix[-2, -3] = b3 / 2
        matrix[-2, -2] = b4 / 2
        matrix[-2, -1] = b2 / 4

        matrix[-1, -2] = b3
        matrix[-1, -1] = b1 / 2

        return csr_matrix(
            matrix
        )


    def make_b3():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = b5
        matrix[0, 1] = b3

        matrix[1, 0] = b3 / 2
        matrix[1, 1] = -8
        matrix[1, 2] = 2

        for i in range(
            2,
            size - 2,
        ):

            matrix[i, i - 1] = 2
            matrix[i, i] = -8
            matrix[i, i + 1] = 2

        matrix[-2, -3] = 2
        matrix[-2, -2] = -8
        matrix[-2, -1] = b3 / 2

        matrix[-1, -2] = b3
        matrix[-1, -1] = b5

        return csr_matrix(
            matrix
        )


    def make_d1():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = -b1
        matrix[0, 1] = b1
        matrix[0, 2] = b6

        matrix[1, 0] = b1 / 2
        matrix[1, 1] = b7
        matrix[1, 2] = b5
        matrix[1, 3] = b6 / 2

        for i in range(
            2,
            size - 2,
        ):

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

        return csr_matrix(
            matrix
        )


    def make_d2():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = b7
        matrix[0, 1] = b4
        matrix[0, 2] = 2

        matrix[1, 0] = b4 / 2
        matrix[1, 1] = 18
        matrix[1, 2] = -8
        matrix[1, 3] = 1

        for i in range(
            2,
            size - 2,
        ):

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

        return csr_matrix(
            matrix
        )


    def make_d3():

        matrix = np.zeros(
            (size, size)
        )

        matrix[0, 0] = b8
        matrix[0, 1] = b4
        matrix[0, 2] = 2

        matrix[1, 0] = b4 / 2
        matrix[1, 1] = 19
        matrix[1, 2] = -8
        matrix[1, 3] = 1

        for i in range(
            2,
            size - 2,
        ):

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

        return csr_matrix(
            matrix
        )


    # Diagonal matrices
    I1 = np.zeros(
        (size, size)
    )

    I1[0, 0] = b6
    I1[-1, -1] = b6

    for i in range(
        1,
        size - 1,
    ):

        I1[i, i] = 2


    I1 = csr_matrix(
        I1
    )

    I2 = I1 / 2


    B1 = make_b1()
    B2 = make_b2()
    B3 = make_b3()

    D1 = make_d1()
    D2 = make_d2()
    D3 = make_d3()


    blocks = []

    for row in range(size):

        blocks.append(
            [None] * size
        )


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
    for row in range(
        2,
        size - 2,
    ):

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


    return bmat(
        blocks,
        format="csr",
    )


# -------------------------
# Build the plate matrix
# -------------------------

print(
    "Building free-edge matrix..."
)

M = build_free_edge_matrix()

print(
    "Matrix finished."
)


# -------------------------
# Find the natural frequencies
# -------------------------

print(
    "Calculating natural frequencies..."
)


eigenvalues = eigs(
    M,
    k=number_of_eigenvalues,
    sigma=1e-8,
    which="LM",
    return_eigenvectors=False,
)


frequencies = []

all_positive_frequencies = []


for eigenvalue in eigenvalues:

    # The physical eigenvalues should
    # be essentially real.
    if abs(
        eigenvalue.imag
    ) > 1e-7:

        continue


    eigenvalue = (
        eigenvalue.real
    )


    # Ignore rigid-body zero modes
    # and tiny numerical negatives.
    if eigenvalue <= 1e-12:

        continue


    # lambda = h^4 * omega^2
    omega = (
        np.sqrt(
            eigenvalue
        )
        / h**2
    )


    # Convert to physical Hz
    frequency_hz = (
        D
        * omega
        / (
            2
            * np.pi
            * a**2
        )
    )


    all_positive_frequencies.append(
        frequency_hz
    )


    if (
        minimum_frequency
        <= frequency_hz
        <= maximum_frequency
    ):

        frequencies.append(
            frequency_hz
        )


# Sort from low to high
frequencies.sort()


# -------------------------
# Remove repeated frequencies
# -------------------------

# Square plates can have multiple modes
# with essentially the same frequency.
#
# For this center-forced catalog we only
# need to solve once at that frequency.

resonance_frequencies = []


for frequency in frequencies:

    if len(
        resonance_frequencies
    ) == 0:

        resonance_frequencies.append(
            frequency
        )

        continue


    if (
        abs(
            frequency
            - resonance_frequencies[-1]
        )
        > frequency_tolerance
    ):

        resonance_frequencies.append(
            frequency
        )


print()

print(
    f"Predicted resonances from "
    f"{minimum_frequency} to "
    f"{maximum_frequency} Hz"
)

print(
    "----------------------------------------"
)


for number, frequency in enumerate(
    resonance_frequencies,
    start=1,
):

    print(
        f"{number:2d}: "
        f"{frequency:7.1f} Hz"
    )


# Check that enough eigenvalues
# were calculated.
if len(
    all_positive_frequencies
) > 0:

    highest_frequency = max(
        all_positive_frequencies
    )

    print()

    print(
        f"Highest frequency calculated: "
        f"{highest_frequency:.1f} Hz"
    )


    if (
        highest_frequency
        < maximum_frequency
    ):

        print()

        print(
            "Increase number_of_eigenvalues "
            "because the calculation did not "
            "reach 2000 Hz."
        )


# -------------------------
# Prepare center forcing
# -------------------------

plate_operator = (
    M / h**4
)


I = identity(
    number_of_points,
    format="csr",
)


# Discrete point force at center
force = np.zeros(
    number_of_points
)

force[
    center_index
] = (
    1 / h**2
)


# -------------------------
# Calculate forced response
# near every resonance
# -------------------------

results = {}


print()

print(
    "Calculating center-forced "
    "resonance patterns..."
)

print()


for number, resonance_frequency in enumerate(
    resonance_frequencies,
    start=1,
):

    # Solve just slightly above the
    # ideal resonance so the matrix
    # is not exactly singular.
    driving_frequency = (
        resonance_frequency
        + resonance_offset
    )


    print(
        f"{number:2d} / "
        f"{len(resonance_frequencies)}: "
        f"{resonance_frequency:.1f} Hz"
    )


    # Convert physical Hz into the
    # dimensionless angular frequency
    # used by the model.
    omega0 = (
        2
        * np.pi
        * a**2
        * driving_frequency
        / D
    )


    # Harmonic response:
    #
    # (L - omega^2 I) F = Q
    system_matrix = (
        plate_operator
        - omega0**2 * I
    )


    response = spsolve(
        system_matrix,
        force,
    )


    # Normalize because we care about
    # the nodal shape here.
    maximum_response = np.max(
        np.abs(
            response
        )
    )


    if maximum_response != 0:

        response = (
            response
            / maximum_response
        )


    # Turn the flattened vector
    # back into the square plate.
    response_shape = response.reshape(
        (size, size),
        order="F",
    )


    results[
        resonance_frequency
    ] = response_shape


print()

print(
    "Forced response calculation finished."
)


# -------------------------
# Save individual plots
# -------------------------

for number, resonance_frequency in enumerate(
    resonance_frequencies,
    start=1,
):

    response = results[
        resonance_frequency
    ]


    fig, ax = plt.subplots(
        figsize=(6, 6)
    )


    ax.contour(
        X,
        Y,
        response.T,
        levels=[0],
        colors="black",
    )


    ax.set_title(
        f"{resonance_frequency:.1f} Hz"
    )


    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_aspect(
        "equal"
    )


    plt.tight_layout()


    fig.savefig(
        "plots/center_forced_modes/"
        f"mode_{number:02d}_"
        f"{resonance_frequency:.1f}hz.png",
        dpi=300,
    )


    plt.close(
        fig
    )


# -------------------------
# Combined 3 x 3 plots
# -------------------------

modes_per_figure = 9

columns = 3
rows = 3


number_of_modes = len(
    resonance_frequencies
)


number_of_figures = (
    number_of_modes
    + modes_per_figure
    - 1
) // modes_per_figure


for figure_number in range(
    number_of_figures
):

    start = (
        figure_number
        * modes_per_figure
    )

    end = min(
        start + modes_per_figure,
        number_of_modes,
    )


    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(11, 11),
    )


    axes = np.array(
        axes
    ).reshape(-1)


    for plot_number, mode_index in enumerate(
        range(start, end)
    ):

        ax = axes[
            plot_number
        ]


        resonance_frequency = (
            resonance_frequencies[
                mode_index
            ]
        )


        response = results[
            resonance_frequency
        ]


        ax.contour(
            X,
            Y,
            response.T,
            levels=[0],
            colors="black",
        )


        ax.set_title(
            f"{resonance_frequency:.1f} Hz",
            fontsize=11,
            pad=6,
        )


        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)

        ax.set_aspect(
            "equal"
        )


        ax.set_xticks(
            [-1, 0, 1]
        )

        ax.set_yticks(
            [-1, 0, 1]
        )


    # Hide unused spaces on
    # the final page.
    for plot_number in range(
        end - start,
        len(axes),
    ):

        axes[
            plot_number
        ].axis(
            "off"
        )


    plt.subplots_adjust(
        hspace=0.45,
        wspace=0.35,
    )


    fig.savefig(
        "plots/center_forced_modes/"
        f"center_forced_modes_"
        f"part{figure_number + 1}.png",
        dpi=300,
        bbox_inches="tight",
    )


plt.show()
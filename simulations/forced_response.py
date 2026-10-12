"""Center-forced response and nodal lines at selected experimental frequencies."""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve

from plate_model import build_free_edge_matrix


# -------------------------
# Physical plate parameters
# -------------------------

plate_side_length = 17 * 0.0254      # 17-inch square plate, meters
thickness = (1 / 16) * 0.0254        # 1/16 inch, meters
youngs_modulus = 69e9                # Aluminum, Pa
density = 2700                       # Aluminum, kg/m^3
mu = 0.33                            # Poisson's ratio

a = plate_side_length / 2           # Half the side length

# Convert between physical frequency and the model's dimensionless frequency
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

N = 200                            # Use 40 for a faster test
h = 2 / N
size = N + 1
number_of_points = size**2

x = np.linspace(-1, 1, size)
y = np.linspace(-1, 1, size)
X, Y = np.meshgrid(x, y)

center = N // 2
center_index = center + center * size

output_folder = "plots/v2/forced_response"
os.makedirs(output_folder, exist_ok=True)


# -------------------------
# Solve the center-forced response
# -------------------------

print("Building free-edge matrix...")
M = build_free_edge_matrix(N, mu)
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

        # Undamped steady-state response at the selected frequency:
        # (plate_operator - omega0^2 * identity) response = force
        system = plate_operator - omega0**2 * identity
        response = spsolve(system, force)

        # Convert the vector back to a square plate displacement field.
        u = response.reshape((size, size), order="F")
        results[label] = u
        print(f"  Finished {label} Hz")

        # Individual nodal-line plot.
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

    # Turn off unused panels.
    for index in range(len(labels), len(axes)):
        axes[index].axis("off")

    fig.suptitle(group_name.replace("_", " ").title(), fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(folder, "combined.png"), dpi=300)

    # Keep the combined figures open for comparison.
    print(f"Saved plots in {folder}")

plt.show()
print(f"\nPlots saved in {output_folder}")
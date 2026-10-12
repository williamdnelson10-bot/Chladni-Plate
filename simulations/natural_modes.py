
"""Natural frequencies and nodal lines for a free-edge square aluminum plate."""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigs

from plate_model import build_free_edge_matrix


# -------------------------
# Physical plate parameters
# -------------------------

plate_side_length = 17 * 0.0254    # 17 inches, meters
thickness = (1 / 16) * 0.0254     # 1/16 inch, meters

youngs_modulus = 69e9             # Aluminum, Pa
density = 2700                     # kg/m^3
mu = 0.33                         # Poisson's ratio

a = plate_side_length / 2        # Half of the plate side length
D = np.sqrt(
    youngs_modulus * thickness**2
    / (12 * density * (1 - mu**2))
)


# -------------------------
# Frequency range and grid
# -------------------------

minimum_frequency = 50
maximum_frequency = 1000
number_of_eigenvalues = 100
frequency_tolerance = 0.1         # Hz, groups nearly equal frequencies

N = 200                            # Use 40 for a faster test
h = 2 / N
size = N + 1

x = np.linspace(-1, 1, size)
y = np.linspace(-1, 1, size)
X, Y = np.meshgrid(x, y)

center = N // 2
center_index = center + center * size

project_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_folder = os.path.join(project_folder, "plots", "natural_modes")
os.makedirs(output_folder, exist_ok=True)


# -------------------------
# Calculate natural modes
# -------------------------

print("Building free-edge matrix...")
M = build_free_edge_matrix(N, mu)
print("Matrix finished.")

print("Calculating natural modes...")
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

    # Ignore complex values and the near-zero rigid-body modes.
    if abs(eigenvalue.imag) > 1e-7:
        continue
    eigenvalue = eigenvalue.real
    if eigenvalue <= 1e-12:
        continue

    # Convert the eigenvalue into a physical frequency in Hz.
    omega = np.sqrt(eigenvalue) / h**2
    frequency_hz = D * omega / (2 * np.pi * a**2)
    all_positive_frequencies.append(frequency_hz)

    if not (minimum_frequency <= frequency_hz <= maximum_frequency):
        continue

    # Make the eigenvector real and scale its largest displacement to 1.
    mode_vector = eigenvectors[:, index]
    largest_index = np.argmax(np.abs(mode_vector))
    phase = np.angle(mode_vector[largest_index])
    mode_vector = (mode_vector * np.exp(-1j * phase)).real

    maximum_amplitude = np.max(np.abs(mode_vector))
    if maximum_amplitude != 0:
        mode_vector = mode_vector / maximum_amplitude

    # Amount of relative motion at the plate's center.
    center_response = abs(mode_vector[center_index])

    # Change the one-dimensional eigenvector into a square grid.
    mode_shape = mode_vector.reshape((size, size), order="F")

    modes.append({
        "frequency": frequency_hz,
        "center_response": center_response,
        "shape": mode_shape,
    })

modes.sort(key=lambda mode: mode["frequency"])


# -------------------------
# Group repeated frequencies
# -------------------------

frequency_groups = []

for mode in modes:
    if not frequency_groups:
        frequency_groups.append([mode])
        continue

    previous_frequency = frequency_groups[-1][0]["frequency"]

    if abs(mode["frequency"] - previous_frequency) <= frequency_tolerance:
        frequency_groups[-1].append(mode)
    else:
        frequency_groups.append([mode])

# As in the original program, show one representative per group:
# choose the eigenvector with the most motion at the center.
natural_modes = []

for group in frequency_groups:
    representative = max(group, key=lambda mode: mode["center_response"])
    representative["group_size"] = len(group)
    natural_modes.append(representative)


# -------------------------
# Print natural frequencies
# -------------------------

print(f"\nPredicted natural modes: {minimum_frequency}-{maximum_frequency} Hz")
print("--------------------------------------------------")

for mode_number, mode in enumerate(natural_modes, start=1):
    frequency = mode["frequency"]
    center_percent = 100 * mode["center_response"]
    group_size = mode["group_size"]

    repeated_text = ""
    if group_size > 1:
        repeated_text = f"({group_size} modes share this frequency)"

    print(
        f"{mode_number:2d}: {frequency:7.1f} Hz   "
        f"Center = {center_percent:5.1f}%   {repeated_text}"
    )

if all_positive_frequencies:
    highest_frequency = max(all_positive_frequencies)
    print(f"\nHighest frequency calculated: {highest_frequency:.1f} Hz")

    if highest_frequency < maximum_frequency:
        print("Increase number_of_eigenvalues to cover the full range.")


# -------------------------
# Save individual nodal plots
# -------------------------

for mode_number, mode in enumerate(natural_modes, start=1):
    frequency = mode["frequency"]
    mode_shape = mode["shape"]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contour(X, Y, mode_shape.T, levels=[0], colors="black")
    ax.set_title(f"{frequency:.1f} Hz")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect("equal")
    plt.tight_layout()

    filename = f"mode_{mode_number:02d}_{frequency:.1f}hz.png"
    fig.savefig(os.path.join(output_folder, filename), dpi=300)
    plt.close(fig)


# -------------------------
# Save combined 3 x 3 plots
# -------------------------

modes_per_figure = 9
number_of_figures = (len(natural_modes) + modes_per_figure - 1) // modes_per_figure

for figure_number in range(number_of_figures):
    start = figure_number * modes_per_figure
    end = min(start + modes_per_figure, len(natural_modes))

    fig, axes = plt.subplots(3, 3, figsize=(11, 11))
    axes = axes.reshape(-1)

    for plot_number, mode_index in enumerate(range(start, end)):
        ax = axes[plot_number]
        mode = natural_modes[mode_index]

        ax.contour(X, Y, mode["shape"].T, levels=[0], colors="black")
        ax.set_title(f"{mode['frequency']:.1f} Hz", fontsize=11, pad=6)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect("equal")
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])

    # Turn off blank panels on the last page.
    for plot_number in range(end - start, len(axes)):
        axes[plot_number].axis("off")

    fig.suptitle(
        f"Predicted Natural Modes: {minimum_frequency}-{maximum_frequency} Hz",
        fontsize=16,
    )
    plt.subplots_adjust(top=0.92, hspace=0.45, wspace=0.35)

    filename = (
        f"natural_modes_{minimum_frequency}_{maximum_frequency}hz_"
        f"part{figure_number + 1}.png"
    )
    fig.savefig(
        os.path.join(output_folder, filename), dpi=300, bbox_inches="tight"
    )
    # Leave the combined figures open, as in the original program.

plt.show()
print(f"\nPlots saved in {output_folder}")

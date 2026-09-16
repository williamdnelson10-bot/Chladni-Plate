import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Number of intervals along each side of the plate
N = 200

# Distance between neighboring grid points
h = 2 / N

# Time step
dt = 0.000025

# Final simulation time
final_time = 0.025

# Number of time steps needed
number_of_steps = round(final_time / dt)

# Coordinates of all grid points from -1 to 1
x = np.linspace(-1, 1, N + 1)
y = np.linspace(-1, 1, N + 1)


def calculate_w(u_values):
    # Create an array to store w = laplacian(u)
    w_values = np.zeros((N + 1, N + 1))

    # Interior grid points
    for i in range(1, N):
        for j in range(1, N):
            w_values[i, j] = (
                u_values[i + 1, j]
                + u_values[i - 1, j]
                + u_values[i, j + 1]
                + u_values[i, j - 1]
                - 4 * u_values[i, j]
            ) / h**2

    # Left and right edges
    for j in range(1, N):
        w_values[0, j] = (
            2 * u_values[1, j]
            + u_values[0, j + 1]
            + u_values[0, j - 1]
            - 4 * u_values[0, j]
        ) / h**2

        w_values[N, j] = (
            2 * u_values[N - 1, j]
            + u_values[N, j + 1]
            + u_values[N, j - 1]
            - 4 * u_values[N, j]
        ) / h**2

    # Bottom and top edges
    for i in range(1, N):
        w_values[i, 0] = (
            u_values[i + 1, 0]
            + u_values[i - 1, 0]
            + 2 * u_values[i, 1]
            - 4 * u_values[i, 0]
        ) / h**2

        w_values[i, N] = (
            u_values[i + 1, N]
            + u_values[i - 1, N]
            + 2 * u_values[i, N - 1]
            - 4 * u_values[i, N]
        ) / h**2

    # Four corners
    w_values[0, 0] = (
        2 * u_values[1, 0]
        + 2 * u_values[0, 1]
        - 4 * u_values[0, 0]
    ) / h**2

    w_values[0, N] = (
        2 * u_values[1, N]
        + 2 * u_values[0, N - 1]
        - 4 * u_values[0, N]
    ) / h**2

    w_values[N, 0] = (
        2 * u_values[N - 1, 0]
        + 2 * u_values[N, 1]
        - 4 * u_values[N, 0]
    ) / h**2

    w_values[N, N] = (
        2 * u_values[N - 1, N]
        + 2 * u_values[N, N - 1]
        - 4 * u_values[N, N]
    ) / h**2

    return w_values


def calculate_acceleration(w_values):
    # Create an array to store u_tt
    acceleration_values = np.zeros((N + 1, N + 1))

    # Interior grid points
    for i in range(1, N):
        for j in range(1, N):
            acceleration_values[i, j] = -(
                w_values[i + 1, j]
                + w_values[i - 1, j]
                + w_values[i, j + 1]
                + w_values[i, j - 1]
                - 4 * w_values[i, j]
            ) / h**2

    # Left and right edges
    for j in range(1, N):
        acceleration_values[0, j] = -(
            2 * w_values[1, j]
            + w_values[0, j + 1]
            + w_values[0, j - 1]
            - 4 * w_values[0, j]
        ) / h**2

        acceleration_values[N, j] = -(
            2 * w_values[N - 1, j]
            + w_values[N, j + 1]
            + w_values[N, j - 1]
            - 4 * w_values[N, j]
        ) / h**2

    # Bottom and top edges
    for i in range(1, N):
        acceleration_values[i, 0] = -(
            w_values[i + 1, 0]
            + w_values[i - 1, 0]
            + 2 * w_values[i, 1]
            - 4 * w_values[i, 0]
        ) / h**2

        acceleration_values[i, N] = -(
            w_values[i + 1, N]
            + w_values[i - 1, N]
            + 2 * w_values[i, N - 1]
            - 4 * w_values[i, N]
        ) / h**2

    # Four corners
    acceleration_values[0, 0] = -(
        2 * w_values[1, 0]
        + 2 * w_values[0, 1]
        - 4 * w_values[0, 0]
    ) / h**2

    acceleration_values[0, N] = -(
        2 * w_values[1, N]
        + 2 * w_values[0, N - 1]
        - 4 * w_values[0, N]
    ) / h**2

    acceleration_values[N, 0] = -(
        2 * w_values[N - 1, 0]
        + 2 * w_values[N, 1]
        - 4 * w_values[N, 0]
    ) / h**2

    acceleration_values[N, N] = -(
        2 * w_values[N - 1, N]
        + 2 * w_values[N, N - 1]
        - 4 * w_values[N, N]
    ) / h**2

    return acceleration_values


# Displacement of the plate at t = 0
u = np.zeros((N + 1, N + 1))

# Calculate the initial displacement at every grid point
for i in range(N + 1):
    for j in range(N + 1):
        u[i, j] = (
            np.cos(np.pi * x[i]) * np.cos(2 * np.pi * y[j])
            + np.cos(2 * np.pi * x[i]) * np.cos(np.pi * y[j])
        )


# Calculate the initial acceleration
w = calculate_w(u)
acceleration = calculate_acceleration(w)

# Save the displacement at t = 0
u_previous = u.copy()

# Calculate the displacement after the first time step
u_current = u + 0.5 * dt**2 * acceleration


# Continue advancing the plate through time
for step in range(1, number_of_steps):
    # Calculate w and acceleration using the current displacement
    w = calculate_w(u_current)
    acceleration = calculate_acceleration(w)

    # Calculate the displacement at the next time step
    u_next = (
        2 * u_current
        - u_previous
        + dt**2 * acceleration
    )

    # Move the stored time levels forward
    u_previous = u_current
    u_current = u_next


# Calculate the exact solution at the final time
u_exact = np.zeros((N + 1, N + 1))

for i in range(N + 1):
    for j in range(N + 1):
        u_exact[i, j] = (
            np.cos(np.pi * x[i]) * np.cos(2 * np.pi * y[j])
            + np.cos(2 * np.pi * x[i]) * np.cos(np.pi * y[j])
        ) * np.cos(5 * np.pi**2 * final_time)


# Find the absolute error
error = np.abs(u_current - u_exact)

print("Maximum error:", np.max(error))


# Create x-y grids for the contour plots
X, Y = np.meshgrid(x, y)

# Create one figure with three plots
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))


# Numerical nodal lines
axes[0].contour(
    X,
    Y,
    u_current.T,
    levels=[0],
    colors="black"
)

axes[0].set_title("Numerical Nodal Lines")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
axes[0].set_aspect("equal")


# Exact nodal lines
axes[1].contour(
    X,
    Y,
    u_exact.T,
    levels=[0],
    colors="black"
)

axes[1].set_title("Exact Nodal Lines")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
axes[1].set_aspect("equal")


# Overlay numerical and exact nodal lines
axes[2].contour(
    X,
    Y,
    u_current.T,
    levels=[0],
    colors="black"
)

axes[2].contour(
    X,
    Y,
    u_exact.T,
    levels=[0],
    colors="red",
    linestyles="dashed"
)

axes[2].set_title("Nodal Line Comparison")
axes[2].set_xlabel("x")
axes[2].set_ylabel("y")
axes[2].set_aspect("equal")


# Legend for the comparison plot
legend_lines = [
    Line2D([0], [0], color="black", label="Numerical"),
    Line2D([0], [0], color="red", linestyle="--", label="Exact")
]

axes[2].legend(handles=legend_lines)


# Title for the whole figure
fig.suptitle(f"Validation Case at t = {final_time}")

plt.tight_layout()
plt.show()
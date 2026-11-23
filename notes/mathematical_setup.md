# Basics
We model a thin square plate with free edges. The square geometry gives useful symmetry, but the free-edge conditions are difficult to solve analytically, so we use a numerical model. A simpler problem with different boundary conditions provides an exact solution for an initial check.
We want to predict the motion of this plate vibrating at a specific frequency, as we drive a force through the center of it.
At resonance, you get the strongest response and most distinct visual pattern visible in nodal lines.
    These nodal lines are where the sand we pour onto the plate will collect, because the anti-nodal regions will vibrate it towards these lines, following whatever natural frequency pattern occurs at that particular mode.
The center driving force creates waves which travel outwards towards the edge then reflect back to interfere with each other, to create standing waves.

We will use the elasticity model, which relates to the deformation of an elastic (the plate in this case) due to stress, strain, and displacement distribution resulting from external forces (driving force, boundary/support forces, etc.)
    So rather than treating the plate like an ideal stretched membrane, we factor in its elastic bending stiffness.

# Research Background
Linqi Shao wrote a research paper in 2018 first modelling the plate as an ideal stretched membrane, using the 2D wave equation. But the nodal patterns didn't line up with the measured patterns. Then they used the model that we'll use for the actual case, which was a fourth order elasticity equation. They got more accurate results, but they had issues dealing with the boundary conditions.
D'Alessio later wrote a research paper published in 2021 using the same equation, but had a better technique for dealing with the troublesome boundary conditions, and got more accurate results. So this paper is what we'll be modelling our approach after.

# Simpler Case
D'Alessio first simulates a simpler case with a known solution, which can be used to compare to the proper initial and boundary conditions, where there is no known solution. In this simpler case he still uses the elasticity equation, but with no force term $Q$:

$$
\frac{\partial^2 u}{\partial t^2}+D^2\nabla^4u=0
$$

$u(x,y,t)$ is the transverse displacement of the plate, $x$ and $y$ are the spatial coordinates and $t$ for time, $D^2$ has to do with the flexural rigidity of the plate, and $\nabla^4u$ is the biharmonic operator which describes the plate's bending response, as a fourth order spatial derivative.

The boundary conditions and initial conditions are chosen such that the equation has an exact solution.

Boundary conditions:

$$
\frac{\partial u}{\partial x}=0,\qquad \frac{\partial w}{\partial x}=0\qquad \text{along }x=\pm a
$$

$$
\frac{\partial u}{\partial y}=0,\qquad \frac{\partial w}{\partial y}=0\qquad \text{along }y=\pm a
$$

$$
w=\nabla^2u=\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}
$$

$w$ is a variable we use temporarily to simplify the simulation process, defined as the Laplacian of the displacement.

However, it's easier to use these equations and initial/boundary conditions, as well as the versions we'll use later in a non-dimensional form:

Equation:

$$
\frac{\partial^2 u}{\partial t^2}+\nabla^4 u=0
$$

Boundary conditions:

$$
\frac{\partial u}{\partial x}=0,\qquad \frac{\partial w}{\partial x}=0\qquad \text{along }x=\pm1
$$

$$
\frac{\partial u}{\partial y}=0,\qquad \frac{\partial w}{\partial y}=0\qquad \text{along }y=\pm1
$$

Initial Conditions:

$$
U_0(x,y)=\cos(\pi x)\cos(2\pi y)+\cos(2\pi x)\cos(\pi y)
$$

$$
\dot{U}_0=0
$$

These dimensionless forms are facilitated by the relations:

$$
x=ax^\ast,\qquad y=ay^\ast,\qquad t=\frac{a^2}{D}t^\ast,\qquad u=Hu^\ast,\qquad w=\frac{H}{a^2}w^\ast,\qquad U_0=HU_0^\ast
$$

With values:

$$
x^\ast=\frac{x}{a},\qquad y^\ast=\frac{y}{a},\qquad t^\ast=\frac{D}{a^2}t,\qquad u^\ast=\frac{u}{H},\qquad w^\ast=\nabla^{\ast 2}u^\ast,\qquad \dot{U}_0=\frac{HD}{a^2}\dot{U}_0^\ast
$$

# Actual Case
In the case which actually matches the conditions of the experiment, we now have an oscillating force at the center. We have the real boundary conditions, which also involve the corners, as well as initial conditions which aren't artificially set at some displacement, there is no initial displacement or initial velocity.

Equation:

$$
\frac{\partial^2u}{\partial t^2}+\nabla^4u=Q(x,y,t)
$$

$$
Q(x,y,t)=\delta(x)\delta(y)\cos(\omega_0t)
$$

(Q(x,y,t)) represents the oscillating force applied at the center of the plate, modeled using Dirac delta functions. Because we're using the dimensionless equation, the angular driving frequency must also be dimensionless. We convert the physical driving frequency (f_0), measured in Hz, using

$$
\omega_0^\ast=\frac{2\pi a^2}{D}f_0
$$

where (a) is half the plate's side length. In the dimensionless equations, we drop the asterisk for simplicity.

Boundary Conditions:

$$
\frac{\partial^2u}{\partial x^2}+\mu\frac{\partial^2u}{\partial y^2}=0\qquad \text{along }x=\pm1
$$

$$
\frac{\partial^3u}{\partial x^3}+(2-\mu)\frac{\partial^3u}{\partial x\partial y^2}=0\qquad \text{along }x=\pm1
$$

$$
\frac{\partial^2u}{\partial y^2}+\mu\frac{\partial^2u}{\partial x^2}=0\qquad \text{along }y=\pm1
$$

$$
\frac{\partial^3u}{\partial y^3}+(2-\mu)\frac{\partial^3u}{\partial x^2\partial y}=0\qquad \text{along }y=\pm1
$$

$$
\frac{\partial^2u}{\partial x\partial y}=0\qquad \text{at }(x,y)=(\pm1,\pm1)
$$

Initial Conditions:

$$
u(x,y,0)=0
$$

$$
\frac{\partial u}{\partial t}(x,y,0)=0
$$

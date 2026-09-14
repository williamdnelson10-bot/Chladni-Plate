# Basics

Our goal is to simulate the motion of a thin square plate with free edges, the free edges and square shape make the math much simpler, giving us a clean analytic solution we can compare our results to.
We want to predict the motion of this plate vibrating at a specific frequency, as we drive a force through the center of it.
At resonance you get the strongest response and most distinct visual pattern in nodal lines.
    These nodal lines are where the sand we pour onto the plate will collect, because the anti-nodal regions will vibrate it towards these lines, following whatever natural frequency pattern occurs at that particular mode.
The center driving force creates waves which travel outwards towards the edge then reflect back to interefere with each other, to create standing waves.

We will use the elasticity model, which relates to the deformation of an elastic (the plate in this case) due to stress, strain, and displacement distribution resulting from external forces (driving force, boundry/support forces, etc.)
    So rather than treating the plate like an ideal stretched membrane, we factor in its elastic bending stiffness.

# Research Background

Linqu Shao wrote a research paper in 2018 first modelling the plate using the 2D wave equation, but the nodal patterns didn't line up with the measured patterns. Then they used the model that we'll use for the actual case, which was a fourth order elasticity equation. They got more accurate results, but they had issues dealing with the boundry conditions.
D'Alessio later wrote a research paper published in 2021 using the same equation, but had a better technique for dealing with the troublesome boundry conditions, and got better results. So this paper is what we'll be modelling our approach after.

# Simpler Case
D'Alessio first simulates a simpler case with a known solution, which can be used to compare to the proper initial and boundry conditions, where there is no known solution. In this simpler case he still uses the elasticity equation, but with no force term Q:

$$
\frac{\partial^2 u}{\partial t^2}+D^2\nabla^4u=0
$$

But with boundry conditions and initial conditions chosen such that the equation has an exact solution.
Boundry conditions:

$$
\frac{\partial u}{\partial x}=0,\qquad \frac{\partial w}{\partial x}=0\qquad \text{along }x=\pm a
$$

$$
\frac{\partial u}{\partial y}=0,\qquad \frac{\partial w}{\partial y}=0\qquad \text{along }y=\pm a
$$

$$
w=\nabla^2u=\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}
$$

However, it's easier to use these equation and boundry conditions, as well as the boundry conditions we'll use later in a non-dimensional form:

Equation:

$$
\frac{\partial^2 u}{\partial t^2}+\nabla^4 u=0
$$

Boundry conditions:

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
\dot U_0=0
$$

These dimensionless forms are facillitated by these relations:

$$
x=ax^\ast,\qquad y=ay^\ast,\qquad t=\frac{a^2}{D}t^\ast,\qquad u=Hu^\ast,\qquad w=\frac{H}{a^2}w^\ast,\qquad U_0=HU_0^\ast
$$

With values:

$$
x^\ast=\frac{x}{a},\qquad y^\ast=\frac{y}{a},\qquad t^\ast=\frac{D}{a^2}t,\qquad u^\ast=\frac{u}{H},\qquad w^\ast=\nabla^{\ast 2}u^\ast,\qquad \dot U_0=\frac{HD}{a^2}\dot U_0^\ast
$$

# Actual Case
In the case which actually matches the conditions of the experiment, we now have an oscillating force at the center, the real plate conditions, which also involve the corners, as well as initial conditions which aren't artifically set at some displacement, there is no initial displacement or initial velocity.

Equation:

$$
\frac{\partial^2u}{\partial t^2}+\nabla^4u=Q(x,y,t)
$$

$$
Q(x,y,t)=\delta(x)\delta(y)\cos(\omega_0t)
$$

Boundry Conditions:

$$
\frac{\partial^2u}{\partial x^2}+\nu\frac{\partial^2u}{\partial y^2}=0\qquad \text{along }x=\pm1
$$

$$
\frac{\partial^3u}{\partial x^3}+(2-\nu)\frac{\partial^3u}{\partial x\partial y^2}=0\qquad \text{along }x=\pm1
$$

$$
\frac{\partial^2u}{\partial y^2}+\nu\frac{\partial^2u}{\partial x^2}=0\qquad \text{along }y=\pm1
$$

$$
\frac{\partial^3u}{\partial y^3}+(2-\nu)\frac{\partial^3u}{\partial x^2\partial y}=0\qquad \text{along }y=\pm1
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

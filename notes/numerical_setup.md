# Intro
We want to simulate the behavior of the thin plate. The method used in D'Alessio's paper is called the finite volume difference. We divide the plate into a grid of points, and model each point as having some displacement u which extends halfway towards the next closest point in every spatial direction. However, all of this spatial information, as well as the acceleration information we will calculate, are all contained at each grid point. The method will be explained in the more realstic experimental case, but the idea is similar for the simpler case, just without a $Q$ term, and with different initial conditions and boundries.

# Finite Volume Method
You divide the plate up into a series of equally spaced points along the $x$ and $y$ directions. The initial displacement is zero, and the oscillations are caused by the center driven force. We turn on the force right at $t = 0$, such that the plate hasn't moved yet, but it will move once we travel some time $\Delta t$, a time step which remains constant throughout this simulation/experiment. At that time, we have some predicted displacement at every point, from the acceleration at the center, as well as from the initial conditions. We use neighboring displacements to estimate the acceleration at every point. We use these new accelerations to calculate a new set of diplacements at $2\Delta t$, then we repeat this process.

However, since the original equation is in the fourth order, we need to use a two step process to make these calculations, using the variable $w$:

$$
\text{1)}\quad w=\nabla^2u
$$

$$
\text{2)}\quad \frac{\partial^2u}{\partial t^2}+\nabla^2w=Q(x,y,t)
$$

# Mathematical Setup
The numerical approximation $u_{i,j}^k$ represents $u(x_i,y_j,t_k)$, at mesh points $x_i=-1+ih$, $y_j=-1+jh$, and $t_k=k\Delta t$, where $h=\frac{2}{N}$.  
D'Alessio uses $h=0.01$ and $\Delta t=0.000025$, so $N=200$.  
Here, $i$ and $j$ are the spatial grid indices, $k$ is the time-step index, $h$ is the spacing between neighboring grid points, $\Delta t$ is the time-step size, and $N$ is the number of intervals along each spatial direction.  
The plate itself uses $-1\leq x,y\leq1$ for its bounds, using the non-dimensional version of the equation and initial conditions.

# Validation Case
We use initial conditions:  
$U_0(x,y)=\cos(\pi x)\cos(2\pi y)+\cos(2\pi x)\cos(\pi y)$  
$\dot U_0=0$  
$u_{i,j}^0=U_0(x_i,y_j)$

## Calculating w
To approximate $w$ at each grid point, we use a five point setup, consisting of the particular point itself, as well as the surrounding points in each cardinal direction. At each interior gridpoint, we don't need to worry about boundry conditions yet, and can approximate $w=\nabla^2u$ using  
$w_{i,j}=\frac{u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j}}{h^2}$  

## Boundry values of w
To get $w$ on the edges, we need to use the boundry conditions:  

$\frac{\partial u}{\partial x}=0$ along $x=\pm1$

$\frac{\partial u}{\partial y}=0$ along $y=\pm1$

We need to use a ghost point just outside the plane to maintain the proper boundry conditions. For example, along the left edge we get:

$u_{-1,j}=u_{1,j}$

$w_{0,j}=\frac{2u_{1,j}+u_{0,j+1}+u_{0,j-1}-4u_{0,j}}{h^2}$

## Corner values of w
Now we're mising neighbors in two directions instead of one, so we apply the ghost point idea in both directions:  

$w_{0,0}=\frac{2u_{1,0}+2u_{0,1}-4u_{0,0}}{h^2}$

This is the bottom left corner (0,0), and it's the same idea for the other corners.

## Acceleration
Now we use the 2nd part of the two step method:

$\frac{\partial^2u}{\partial t^2}+\nabla^2w=0$

There's no $Q$ term, so we rearrange to get:

$u_{tt}=-\nabla^2w$

Now at an interior point:

$u_{tt\,i,j}=-\frac{w_{i+1,j}+w_{i-1,j}+w_{i,j+1}+w_{i,j-1}-4w_{i,j}}{h^2}$

## Acceleration: edges and corners
Along the left and right edges we use boundry condition:

$\frac{\partial w}{\partial x}=0$

Along the bottom and top edges we use:

$\frac{\partial w}{\partial y}=0$

Same idea for ghost points as for $u$, at the left edge:

$w_{-1,j}=w_{1,j}$

With corresponding acceleration:

$u_{tt\,0,j}=-\frac{2w_{1,j}+w_{0,j+1}+w_{0,j-1}-4w_{0,j}}{h^2}$

So at the bottom left corner for example:

$u_{tt\,0,0}=-\frac{2w_{1,0}+2w_{0,1}-4w_{0,0}}{h^2}$

## First time step
We need the displacement at the first new time increment:

$u^1_{i,j}=u^0_{i,j}+\Delta t\,\dot u^0_{i,j}+\frac{1}{2}(\Delta t)^2u_{tt\,i,j}^0$

With the initial velocity of zero, this simplifies to:

$u^1_{i,j}=u^0_{i,j}+\frac{1}{2}(\Delta t)^2u_{tt\,i,j}^0$

## Full time step equation
The second time derivative is approximated by:

$u_{tt\,i,j}^k\approx\frac{u_{i,j}^{k+1}-2u_{i,j}^k+u_{i,j}^{k-1}}{(\Delta t)^2}$

Solving for the next displacement:

$u_{i,j}^{k+1}=2u_{i,j}^k-u_{i,j}^{k-1}+(\Delta t)^2u_{tt\,i,j}^k$

Process for every time step:

$u^k\rightarrow w^k\rightarrow u_{tt}^k\rightarrow u^{k+1}$

## Time Loop

Now that we have both $u^0$ and $u^1$, we can advance the simulation through the remaining time steps.

At each time $t_k$, we first calculate $w^k$ from the current displacement $u^k$, then calculate the acceleration $u_{tt}^k$ from $w^k$.

The next displacement is then calculated using:

$u_{i,j}^{k+1}=2u_{i,j}^k-u_{i,j}^{k-1}+(\Delta t)^2u_{tt\,i,j}^k$.

After calculating $u^{k+1}$, the stored time levels are updated so that the newly calculated solution becomes the current solution, while the current solution becomes the previous solution:

$u^{k+1}\rightarrow u^k,\qquad u^k\rightarrow u^{k-1}$

This process repeats as $u^k\rightarrow w^k\rightarrow u_{tt}^k\rightarrow u^{k+1}$ until the chosen final simulation time is reached.

For this test, $t_{\mathrm{final}}=0.025$. With $\Delta t=0.000025$, this requires 1000 time steps.

## Exact Solution Check

For the chosen initial displacement, the $(1,2)$ and $(2,1)$ modes have the same angular frequency, $\omega_{12}=\omega_{21}=\pi^2(1^2+2^2)=5\pi^2$.

The exact solution is $u_{\mathrm{exact}}(x,y,t)=[\cos(\pi x)\cos(2\pi y)+\cos(2\pi x)\cos(\pi y)]\cos(5\pi^2t)$.

The numerical error at each grid point can be calculated with $e_{i,j}=|u_{i,j}-u_{\mathrm{exact}}(x_i,y_j,t)|$.

# Actual Case
In the proper model of the experiment, we now need to factor in the center driving force, the actual boundry conditions, and the actual inital conditions. We know have an initial displacement $u(x,y,0)=0$, An initial velocity $u_t(x,y,0)=0$, an external force at the center $Q(x,y,t)=\delta(x)\delta(y)\cos(\omega_0t)$, and boundry conditions defined in mathematical_setup.md. 

The new governing equation is:  
$u_{tt}+\nabla^4u=Q(x,y,t)$

## Free Edge Matrix
In the validation case, you could just maintain the boundry conditions with ghost points. But now we need to use a matrix M, which represents the gridpoint plate equation, with the proper free edge conditions incorporated. We use the same square grid paramters as before.

Before using the matrix for the center driven version, we want to make sure it works without $Q$. So we use equation $u_{tt}+\nabla^4u=0$. We assume the plate is vibrating as $u(x,y,t)=F(x,y)\cos(\omega t)$. We model this as an eigenvalue problem with $M\mathbf{F}=\lambda\mathbf{F}$. From this, we get eigenvalue equation: $\lambda=h^4\omega^2$. Here, $M$ is the matrix representation of the plate's spatial equation, $\mathbf{F}$ is the vector containing the displacement pattern $F(x,y)$ at all grid points, $\lambda$ is an eigenvalue of $M$, and $\omega$ is the corresponding angular natural frequency.

Using these eigenvalues we're able to obtain the predicted natural frequencies, then compare them to the numbers reported in the D'Alessio paper. We did all of this because when trying to reproduce their patterns, the results were way off, so this was a method that helped us determine what the specific problem in the code was.

## Driven Harmonic Solution
For the actual case, the plate is driven at a chosen angular frequency $\omega_0$. Rather than stepping forward through time, we look for a harmonic response with the same time dependence as the driving force:

$$
u(x,y,t)=F(x,y)\cos(\omega_0t)
$$

Substituting this into the plate equation separates out the time dependence, leaving a spatial equation that can be solved directly for $F$.

The free-edge matrix $M$ is converted into the spatial plate operator $L=\frac{M}{h^4}$, with resulting gridpoint equation: $(L-\omega_0^2I)\mathbf{F}=\mathbf{Q}$.

$L$ is the discretized spatial plate operator, and $I$ is the identity matrix. Solving this system gives the spatial vibration pattern of the plate at the chosen driving frequency.

## Center Force and Plots
The center force is scaled by $\frac{1}{h^2}$ to represent the point force on the discrete grid. At each chosen frequency, the code builds the matrix:

$$
A=L-\omega_0^2I
$$

and then solves

$$
A\mathbf{F}=\mathbf{Q}
$$

We repeat this process for every frequency we're interested in.

Once the spatial response of $F$ has been solved for a given frequency, the nodel lines are where $F(x,y)=0$. The code creates a contour plot wherever this condition is satisfied, and this is repeated for each relevant frequency. 
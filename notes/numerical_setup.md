# Intro
We want to simulate the behavior of the thin plate. The method used in D'Alessio's paper is called the finite volume difference. We divide the plate into a grid of points, and model each point as having some displacement u which extends halfway towards the next closest point in every spatial direction. However, all of this spatial information, as well as the acceleration information we will calculate, are all contained at each grid point. The method will be explained in the more realstic experimental case, but the idea is similar for the simpler case, just without a Q term, and with different initial conditions and boundries.

# Finite Volume Method
You divide the plate up into a series of equally spaced points along the x and y directions. The initial displacement is zero, and the oscillations are caused by the center driven force. We turn on the force right at t = 0, such that the plate hasn't moved yet, but it will move once we travel some time delta t, a time step which remains constant throughout this simulation/experiment. At that time, we have some predicted displacement at every point, from the acceleration at the center, as well as from the initial conditions. We use neighboring displacements to estimate the acceleration at every point. We use these new accelerations to calculate a new set of diplacements at 2delta t, then we repeat.

However, since the original equation is in the fourth order, we need to use a two step process to make these calculations, using the variable w:
$$
w=\nabla^2u
$$

$$
\frac{\partial^2u}{\partial t^2}+\nabla^2w=Q(x,y,t)
$$

# Mathematical Setup
The numerical approximation $u_{i,j}^k$ represents $u(x_i,y_j,t_k)$, at mesh points $x_i=-1+ih$, $y_j=-1+jh$, and $t_k=k\Delta t$, where $h=\frac{2}{N}$.  
D'Alessio uses $h=0.01$ and $\Delta t=0.000025$, so $N=200$.  
The plate itself uses $-1\leq x,y\leq1$ for its bounds, using the non-dimensional version of the equation and initial conditions.

# Validation Case
We use initial conditions:  
$U_0(x,y)=\cos(\pi x)\cos(2\pi y)+\cos(2\pi x)\cos(\pi y)$  
$\dot U_0=0$  
$u_{i,j}^0=U_0(x_i,y_j)$

## Calculating w
At each interior gridpoint, we approximate $w=\nabla^2u$ using  
$w_{i,j}=\frac{u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j}}{h^2}$  
We can do this because every interior point has neighbors on all sides.

## Boundry values of w
To get w on the edges, we need to use the boundry conditions:  

$\frac{\partial u}{\partial x}=0$ along $x=\pm1$

$\frac{\partial u}{\partial y}=0$ along $y=\pm1$

We need to use a ghost point just outside the plane to maintain the proper boundry conditions:

$u_{-1,j}=u_{1,j}$

$w_{0,j}=\frac{2u_{1,j}+u_{0,j+1}+u_{0,j-1}-4u_{0,j}}{h^2}$

## Corner values of w
Now we're mising neighbors in two directions instead of one, so we apply the ghost point idea in both directions:  

$w_{0,0}=\frac{2u_{1,0}+2u_{0,1}-4u_{0,0}}{h^2}$

This is the bottom left corner (0,0), and it's the same idea for the other corners.

## Acceleration
Now we use the 2nd part of the two step method:

$\frac{\partial^2u}{\partial t^2}+\nabla^2w=0$

We can ignore Q, and rearrange to get:

$u_{tt}=-\nabla^2w$

Now at an interior point:

$u_{tt\,i,j}=-\frac{w_{i+1,j}+w_{i-1,j}+w_{i,j+1}+w_{i,j-1}-4w_{i,j}}{h^2}$

## Acceleration: edges and corners
Along the left and right edges we use boundry condition:

$\frac{\partial w}{\partial x}=0$

Along the bottom and top edges we use:

$\frac{\partial w}{\partial y}=0$

Same idea for ghost points as for u, at the left edge:

$w_{-1,j}=w_{1,j}$

with corresponding acceleration:

$u_{tt\,0,j}=-\frac{2w_{1,j}+w_{0,j+1}+w_{0,j-1}-4w_{0,j}}{h^2}$

So at the bottom left corner for example:

$u_{tt\,0,j}=-\frac{2w_{1,j}+w_{0,j+1}+w_{0,j-1}-4w_{0,j}}{h^2}$

## First time step
We need the displacement at the first new time increment:

$u^1_{i,j}=u^0_{i,j}+\Delta t\,\dot u^0_{i,j}+\frac{1}{2}(\Delta t)^2u_{tt\,i,j}^0$

With the initial velocity of zero, this simplifies to:

$u^1_{i,j}=u^0_{i,j}+\frac{1}{2}(\Delta t)^2u_{tt\,i,j}^0$

## Full time step equation
The second time derivative is approximated by:

$u_{tt\,i,j}^k\approx\frac{u_{i,j}^{k+1}-2u_{i,j}^k+u_{i,j}^{k-1}}{(\Delta t)^2}$

Solving for the next displacement:

$u_{tt\,i,j}^k\approx\frac{u_{i,j}^{k+1}-2u_{i,j}^k+u_{i,j}^{k-1}}{(\Delta t)^2}$

Process for every time step:

$u^k\rightarrow w^k\rightarrow u_{tt}^k\rightarrow u^{k+1}$

## Time Loop

Now that we have both $u^0$ and $u^1$, we can advance the simulation through the remaining time steps.

At each time $t_k$, we first calculate $w^k$ from the current displacement $u^k$, then calculate the acceleration $u_{tt}^k$ from $w^k$.

The next displacement is then calculated using $u_{i,j}^{k+1}=2u_{i,j}^k-u_{i,j}^{k-1}+(\Delta t)^2u_{tt\,i,j}^k$.

After calculating $u^{k+1}$, the stored time levels are shifted forward so that $u^{k-1}\leftarrow u^k$ and $u^k\leftarrow u^{k+1}$.

This process repeats as $u^k\rightarrow w^k\rightarrow u_{tt}^k\rightarrow u^{k+1}$ until the chosen final simulation time is reached.

For this test, $t_{\mathrm{final}}=0.025$. With $\Delta t=0.000025$, this requires 1000 time steps.

## Exact Solution Check

For the chosen initial displacement, the $(1,2)$ and $(2,1)$ modes have the same angular frequency, $\omega_{12}=\omega_{21}=\pi^2(1^2+2^2)=5\pi^2$.

The exact solution is $u_{\mathrm{exact}}(x,y,t)=[\cos(\pi x)\cos(2\pi y)+\cos(2\pi x)\cos(\pi y)]\cos(5\pi^2t)$.

The numerical error at each grid point can be calculated with $e_{i,j}=|u_{i,j}-u_{\mathrm{exact}}(x_i,y_j,t)|$.
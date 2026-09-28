# Physical Frequency Scaling
Recommended to read the mathematical_setup and numerical_setup first

In the D'Alessio paper, they used their experimental parameters to convert the dimensional equations and relevant values into non-dimensional ones. Then they used the natural frequencies associated with these, and that's why the code itself doesn't explicity have all of the relevant quantities mentioned. The quantities which are implicity used in the code are $E,\ \rho,\ l,\ a$. So we need to use the values which pertain to our own setup to determine what frequencies we should use.

## Plate Parameters

The plate parameter $D$ is determined from the material properties and plate thickness:

$$
D=\sqrt{\frac{El^2}{12\rho(1-\mu^2)}}
$$

The material properties for our plate are those associated with aluminum, $E$, $\rho$, and $\mu$. Whereas $l$ is measured directly from the thickness of our plate, and this part is different from what D'Alessio had.

The mathematical model describes the plate over the bounds $-a\leq x,y\leq a$, where $a$ is half of the full side length $L$:

$$
a=\frac{L}{2}
$$

After converting to the dimensionless form, these bounds become $-1\leq x,y\leq1$. So our bounds end up being the same, but the difference in $a$ values is relevant in equation:

$$
f=\frac{D\omega^\ast}{2\pi a^2}
$$

The free edge matrix gives us eigenvalues $\lambda$, which correspond to the different natural modes of the plate. The dimensionless angular frequency is related to the eigenvalue by

$$
\lambda=h^4(\omega^\ast)^2,
\qquad \text{and} \qquad
\omega^\ast=\frac{\sqrt{\lambda}}{h^2}.
$$

## Converting to Physical Frequencies

To get the frequency that we would actually use in the experiment, we convert the dimensionless angular frequency back into a physical frequency.

The physical angular frequency is

$$
\omega=\frac{D}{a^2}\omega^\ast.
$$

Since the experiment uses frequency in Hz rather than angular frequency, we divide by $2\pi$:

$$
f=\frac{\omega_{\mathrm{physical}}}{2\pi}.
$$

Combining these gives

$$
f=\frac{D\omega^\ast}{2\pi a^2}.
$$

So once we have the natural frequencies for the matrix, we can use our plate values to find whatever frequencies we'd want to use in Hz. These frequencies are where we should get the best plate response, giving the most distinct nodal lines. But due to experimental errors, during the actual experiment we'll want to test frequencies around these values to find the strongest responses.
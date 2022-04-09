# BETES

BETES (or **B**aroclinic **E**kman **T**ransport **E**quation **S**olver) is an in-house code developed to simulate wind-driven boundary layer currents by computing the boundary-value problem:
```math
\dfrac{\text{d}}{\text{d}z}\left(-K\dfrac{\text{d}\mathbf{u}}{\text{d}z}\right) = -f\mathbf{u}^\perp - \dfrac{1}{\rho}\nabla p_a - \dfrac{g}{\rho}\int_z^0\nabla\rho\,\text{d}z
```
in a water depth column, $`0\leq z\leq-H`$, where the velocity field is $`\mathbf{u}(z) = (u(z), v(z))^T`$, and $`\mathbf{u}^\perp = (-v(z), u(z))^T`$. The barotropic and baroclinc pressure gradient are assumed $`\nabla p_a = \text{const.}`$ and $`\nabla\rho \propto \exp(z)`$, respectively. Coriolis force, $`f`$, and density, $`\rho`$, are constants.

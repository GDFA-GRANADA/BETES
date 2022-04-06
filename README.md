# BETES

BETES (or **B**aroclinic **E**kman **T**ransport **E**quation **S**olver) is a numerical code to simulate wind-driven boundary layer currents. The BETES code computes the prototype system of ODE:
```math
\dfrac{\text{d}}{\text{d}z}\left(-K\dfrac{\text{d}\mathbf{u}}{\text{d}z}\right) = -f\mathbf{u}^\perp - \dfrac{1}{\rho}\nabla p_a - \dfrac{g}{\rho}\int_z^0\nabla\rho\,\text{d}z + 
```
in a water depth column, $`0\leq z\leq-H`$, where the velocity field is $`\mathbf{u}(z) = (u(z), v(z))^T`$, and $`\mathbf{u}^\perp = (-v(z), u(z))^T`$. The barotropic and baroclinc pressure gradient are assumed $`\nabla p_a`$ and $`\nabla\rho \propto \exp(z)`$, respectively. Coreolis force, $`f`$, and density, $'\rho'$, are constants.    

# BETES

BETES (or **B**aroclinic **E**kman **T**ransport **E**quation **S**olver) is a numerical code to simulate wind-driven boundary layer currents. The BETES code computes the prototype system of ODE:
```math
\mathbf{0} = -f\mathbf{u}^\perp - \dfrac{1}{\rho}\nabla p_a - \dfrac{g}{\rho}\int_z^0\nabla\rho\,\text{d}z + \dfrac{\text{d}}{\text{d}z}\left(K\dfrac{\text{d}\mathbf{u}}{\text{d}z}\right)
```
in a water depth column, $`0\leq z\leq-H`$. 

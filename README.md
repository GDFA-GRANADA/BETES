# Introduction
BETES (or **B**aroclinic **E**kman **T**ransport **E**quation **S**olver) is a wind-driven ocean currents simulator in a surface boundary layer. The BETES code computes the boundary-value problem:
```math
\dfrac{\text{d}}{\text{d}z}\left(-K\dfrac{\text{d}\mathbf{u}}{\text{d}z}\right) = -f\mathbf{u}^\perp - \dfrac{1}{\rho}\nabla p_a - \dfrac{g}{\rho}\int_z^0\nabla\rho\,\text{d}z
```
in a water depth column, $`0\leq z\leq-H`$, where the velocity field is $`\mathbf{u}(z) = (u(z), v(z))^T`$, and $`\mathbf{u}^\perp = (-v(z), u(z))^T`$. The barotropic and baroclinc pressure gradient are assumed $`\nabla p_a = \text{const.}`$ and $`\nabla\rho \propto \exp(z)`$, respectively. Coriolis force, $`f`$, and density, $`\rho`$, are constants.

# Getting started
Several scripts and documents are provided to learn to work with BETES. In folder _src_ we find the `Ekman_model.m` archive which is the main code of BETES. 
On the other hand `adjustpdfpage.m` and `profile_Kz.m` are functions to print the plots in pdf and provide, for the main code, the eddy turbulent viscosity, respectively.
The archive `datos_GoC_medios.mat` contains experimental results taken from Mazagón (Huelva, Spain).

In folder _doc_ we find a pdf to describe roughly the numerical scheme used for BETES.

# Programming language and installation
You need MATLAB (R2020b or later) to run BETES.

# List of Contributors

The list of the people who contributed to BETES are the following:
- **Victor Javier Llorente Lázaro** (Technical University of Madrid)
- **Manuel Diez Minguito** (University of Granada)
- **Enrique Manuel Padilla de la Torre** (University of Granada)


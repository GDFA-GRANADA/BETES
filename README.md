# Introduction
BETES (or **B**aroclinic **E**kman **T**ransport **E**quation **S**olver) is a wind-driven ocean currents simulator in a surface boundary layer. The BETES code computes a steady-state Boussinesq flow with a small Rossby number in a water finite depth column. Coriolis force and Barothopic & Baroclinic pressure gradients are balanced with vertical diffusion momentum. Ocean stratification is considered by a family of shapes for the eddy viscosity profile with a two-layer formulation, and the horizontal density gradient is known and vary in depth. 

# Getting started
Several scripts and documents are provided to learn to work with BETES. In folder _src_ we find the `Ekman_model.m` archive which is the main code of BETES. 
On the other hand `adjustpdfpage.m` and `profile_Kz.m` are functions to print the plots in pdf and provide, for the main code, the eddy turbulent viscosity, respectively.
The archive `datos_GoC_medios.mat` contains experimental results taken from Mazagón (Huelva, Spain).

In folder _doc_ we find a pdf to describe roughly the numerical scheme used for BETES.

# Programming language and installation
You need MATLAB (version R2020b or later) to run BETES.

# List of Contributors
The list of the people who contributed to BETES are the following:
- **Victor Javier Llorente Lázaro** (Technical University of Madrid)
- **Manuel Diez Minguito** (University of Granada)
- **Enrique Manuel Padilla de la Torre** (University of Granada)

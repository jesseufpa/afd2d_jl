# AFD2d_mod.jl Documentation


## Finite-difference modeling for a constant density acoustic wave equation

```math
\frac{1}{c^2(\mathbf{x})} \frac{\partial P}{\partial t} = \nabla^2 P -
\rho \, \dot{q}(t)
\, \delta(\mathbf{x} - \mathbf{x}_{s})
```
where $c(\mathbf{x})$ is the medium propagation velocity, $\rho$ its density,  and
$P(\mathbf{x},t)$ the acoustic wavefield excited by the
volume injection rate $q(t)$ at the source position $\mathbf{x}_{s}$.

#### Authors: Jesse Costa and Joerg Schleicher 

```@meta
CurrentModule = AFD2d_mod
```

```@autodocs
Modules = [AFD2d_mod,SU_mod]
Order = [:function,:type]
``` 


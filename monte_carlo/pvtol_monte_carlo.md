# Monte Carlo simulation and  minimum sample size for a PVTOL 

This example point to the MEOSIA controsl exam in current repo.
The scope is shows how to estimate a more realistic performance metric of a controller PVTOL aircrfact under uncertainty, and how many monte carlo runs are needed to report a test result with **95 % confidence** and an
**error smaller than 0.01**.

Example Script: [pvtol_monte_carlo.py](pvtol_monte_carlo.py)

Requirements: Python 3, `numpy`, optionally `matplotlib` (figure generation).

# PVTOL system

PVTOL model:

```math
\ddot{x} = -\left(\frac{U}{m}\right)\sin(\theta)
```

```math
\ddot{z} = \left(\frac{U}{m}\right)\cos(\theta) - g
```

```math
\ddot{\theta} = \frac{1}{J}\tau
```

For the example lets take the wind force $d = (d_x, d_z)$ is added as an external
disturbance ($+d_x/m$ and $+d_z/m$ in the translational equations).

## Wind disturbance. 

In the real world this kind of model will interact with external forces and unexpected disturbances, the model above is ideal.
To also include environmental uncertainty in simulation, a wind force $d = (d_x, d_z)$ [N] is added. The model equations
come from Newton's second law $m\ddot{\mathbf r} = \sum \mathbf F$, so the wind
is added as one more force and then divided by $m$, like the thrust:

```math
m\ddot{x} = -U\sin(\theta) + d_x \;\Rightarrow\; \ddot{x} = -\frac{U}{m}\sin(\theta) + \frac{d_x}{m}
```

```math
m\ddot{z} = U\cos(\theta) - mg + d_z \;\Rightarrow\; \ddot{z} = \frac{U}{m}\cos(\theta) - g + \frac{d_z}{m}
```
The wind is assumed to act at the center of mass, so it produces no torque and
$\ddot{\theta} = \tau/J$ is unchanged (a wind torque would add a $d_\theta/J$ term).

 
## Parameter values form the exam:

| Symbol | Meaning | Nominal |
|---|---|---|
| $x, z$ | translational positions [m] | state |
| $\theta$ | attitude angle [rad] | state |
| $m$ | vehicle mass [kg] | 1 |
| $J$ | moment of inertia [kg m²] | 1 |
| $g$ | gravity [m/s²] | 9.81 |
| $U$ | total thrust of the motors [N] | control |
| $\tau$ | attitude torque [N m] | control |

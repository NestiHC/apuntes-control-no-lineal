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

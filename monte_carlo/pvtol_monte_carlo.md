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

**Task:** start near the origin and hover at $(x_r, z_r) = (1, 1)$ m, simulated
for $T = 8$ s with $\Delta t = 0.01$ s (semi-implicit Euler).

### Controller

A cascaded PD controller designed with the **nominal** parameters only, so
parameter uncertainty produces real model mismatch:

```math
a_x^{des} = k_{px}(x_r - x) - k_{dx}\dot{x}
```

```math
a_z^{des} = k_{pz}(z_r - z) - k_{dz}\dot{z}
```

```math
\theta^{des} = -\frac{a_x^{des}}{g_0}
```

(since $\ddot{x} \approx -g\theta$ near hover)

```math
U = \frac{m_0 (g_0 + a_z^{des})}{\cos(\theta)}
```

```math
\tau = J_0 \left[ k_{p\theta}(\theta^{des} - \theta) - k_{d\theta}\dot{\theta} \right]
```

The controller uses **noisy measurements** of all states.

## 2. Uncertainty sources (variances)

Each source is modeled as Gaussian $\mathcal N(\mu, \sigma^2)$:

| Group | Source | $\sigma$ | When it is sampled |
|---|---|---|---|
| Parameters | $m$, $J$, $g$ | 0.05 (5 %), 0.1 (10 %), 0.005 | once per run |
| Initial state | $x_0, z_0$ / $\theta_0$ / $\dot x_0, \dot z_0, \dot\theta_0$ | 0.05 / 0.02 / 0.02 | once per run |
| Sensors | position / angle / velocity / angular rate | 0.01 / 0.005 / 0.02 / 0.01 | every step |
| Actuators | $U$ (relative) / $\tau$ (additive) | 2 % / 0.04 N m | every step |
| Wind | constant bias / gust | 0.03 N / 0.05 N | per run / every step |

All values live in the `SIGMA` dictionary of the script and can be changed.

## 3. Output metrics

For each run $i$ the script computes the final position error

```math
e_i = \sqrt{(x_i(T)-x_r)^2 + (z_i(T)-z_r)^2}
```

and two quantities are estimated:

- **Mean error** $\mu_e = E[e]$ (continuous metric).
- **Success probability** $p = P(e < 0.05\text{ m})$ (Bernoulli metric).

## 4. Minimum sample size

For a two-sided 95 % confidence level $z = 1.96$ and a maximum error
(half-width) $E = 0.01$:

**Note on $z$.** In this section $z$ is **not** the altitude of the PVTOL. It is
> the critical value (quantile) of the **standard normal distribution**
> $Z \sim \mathcal N(0, 1)$ that leaves 95 % of the probability between $-z$ and $z$:
>
> ```math
> P(-z \le Z \le z) = 1 - \alpha = 0.95 \;\Rightarrow\; z = z_{1-\alpha/2} = z_{0.975} = 1.96
> ```
>
> It comes from the Central Limit Theorem: for large $n$ the sample mean (and
> the sample proportion) is approximately normal, so its standardized error
> follows $\mathcal N(0, 1)$. Other confidence levels use other quantiles,
> e.g. 90 % → 1.645, 99 % → 2.576. In the script it is `Z_95`.
>
> **From the z-table.** The standard normal table gives the cumulative
> probability $\Phi(z) = P(Z \le z)$. For 95 % two-sided, 2.5 % is left in each
> tail, so look for $\Phi(z) = 0.975$ inside the table: it is in row **1.9**,
> column **0.06**, so $z = 1.9 + 0.06 = 1.96$.
>
> | $z$ | 0.04 | 0.05 | **0.06** | 0.07 |
> |---|---|---|---|---|
> | **1.9** | 0.9738 | 0.9744 | **0.9750** | 0.9756 |
>
> In Python: `scipy.stats.norm.ppf(0.975)` → 1.959964.

**Mean of a continuous variable** (CLT), from the interval $\bar{e} \pm z \sigma / \sqrt{n}$:

```math
n \geq \left( \frac{z \sigma}{E} \right)^2
```

**Proportion** (normal approximation of the binomial), from the interval $\hat{p} \pm z \sqrt{p(1-p)/n}$:

```math
n \geq \frac{z^2 p (1-p)}{E^2}
```

Worst case $p = 0.5$: $n = 1.96^2 \cdot 0.25 / 0.01^2 = 9604$.

Since $\sigma$ and $p$ are unknown, they are first estimated with a **pilot**
Monte Carlo ($n_0 = 500$), then the required $N$ is computed. After the full
run the half-width is recomputed; if it is still above $E$ (pilot estimates are
not exact), extra runs are added.

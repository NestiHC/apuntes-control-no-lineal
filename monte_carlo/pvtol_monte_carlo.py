"""Monte Carlo simulation and minimum sample size for a PVTOL hover controller.

Steps shown by the script:
  1. PVTOL model and nominal parameters.
  2. Uncertainty sources (all variances): parameters, initial conditions,
     sensor noise, actuator noise and wind disturbance.
  3. Nominal simulation (no noise).
  4. Pilot Monte Carlo to estimate sigma (mean metric) and p (success rate).
  5. Minimum sample size for 95 % confidence and error E < 0.01.
  6. Full Monte Carlo with that N.
  7. Verification of the confidence intervals.

PVTOL model (m = 1, J = 1), plus a wind disturbance d = (dx, dz) for the MC study:
  x''  = -(U/m)*sin(th) + dx/m
  z''  =  (U/m)*cos(th) - g + dz/m
  th'' = tau/J
"""
import math

import numpy as np

Z_95 = 1.959964  # two-sided 95 % normal quantile
E_MAX = 0.01     # maximum allowed half-width of the confidence interval

# ---------------------------------------------------------------- Step 1
NOMINAL = {"m": 1.0, "J": 1.0, "g": 9.81}
GAINS = {"kpx": 1.0, "kdx": 2.0, "kpz": 4.0, "kdz": 4.0, "kpt": 25.0, "kdt": 10.0}
TARGET = (1.0, 1.0)          # hover set point (x_r, z_r) [m]
DT, T_FINAL = 0.01, 8.0      # integration step and horizon [s]
SUCCESS_TOL = 0.05           # run is a "success" if final position error < 5 cm

# ---------------------------------------------------------------- Step 2
# Standard deviations; variance = sigma**2.
SIGMA = {
    # plant parameters (true plant differs from controller model)
    "m": 0.05, "J": 0.1, "g": 0.005,
    # initial conditions around the origin
    "x0": 0.05, "z0": 0.05, "th0": 0.02, "vx0": 0.02, "vz0": 0.02, "w0": 0.02,
    # sensor noise (per sample, white)
    "s_pos": 0.01, "s_th": 0.005, "s_vel": 0.02, "s_w": 0.01,
    # actuator noise: U relative, tau additive [N m]
    "a_U": 0.02, "a_tau": 0.04,
    # wind: constant bias per run + white gust per step [N]
    "wind_bias": 0.03, "wind_gust": 0.05,
}


def print_variance_table():
    print("\nSTEP 2 - Uncertainty sources")
    print(f"  {'source':<10}{'sigma':>10}{'variance':>14}")
    for k, s in SIGMA.items():
        print(f"  {k:<10}{s:>10.4f}{s * s:>14.2e}")


def sample_plant(n, rng, scale):
    """Draw n random plants and initial states (scale=0 gives the nominal case)."""
    s = {k: v * scale for k, v in SIGMA.items()}
    p = {k: NOMINAL[k] + s[k] * rng.standard_normal(n) for k in NOMINAL}
    p["state"] = np.vstack([
        s["x0"] * rng.standard_normal(n), s["z0"] * rng.standard_normal(n),
        s["th0"] * rng.standard_normal(n), s["vx0"] * rng.standard_normal(n),
        s["vz0"] * rng.standard_normal(n), s["w0"] * rng.standard_normal(n),
    ])
    p["wind"] = s["wind_bias"] * rng.standard_normal((2, n))
    return p, s


def simulate(n, rng, scale=1.0, record=False):
    """Vectorized simulation of n PVTOL runs; returns final position error per run."""
    p, s = sample_plant(n, rng, scale)
    x, z, th, vx, vz, w = p["state"]
    m0, J0, g0 = NOMINAL["m"], NOMINAL["J"], NOMINAL["g"]
    G = GAINS
    xr, zr = TARGET
    steps = int(T_FINAL / DT)
    hist = np.zeros((steps, 3, n)) if record else None

    for k in range(steps):
        # sensor measurements
        xm = x + s["s_pos"] * rng.standard_normal(n)
        zm = z + s["s_pos"] * rng.standard_normal(n)
        thm = th + s["s_th"] * rng.standard_normal(n)
        vxm = vx + s["s_vel"] * rng.standard_normal(n)
        vzm = vz + s["s_vel"] * rng.standard_normal(n)
        wm = w + s["s_w"] * rng.standard_normal(n)

        # cascaded controller designed with nominal parameters
        ax_des = G["kpx"] * (xr - xm) - G["kdx"] * vxm
        az_des = G["kpz"] * (zr - zm) - G["kdz"] * vzm
        th_des = np.clip(-ax_des / g0, -0.5, 0.5)  # x'' ~ -g*th near hover
        U = m0 * (g0 + az_des) / np.maximum(np.cos(thm), 0.5)
        U = np.clip(U, 0.0, 3.0 * m0 * g0)
        tau = J0 * (G["kpt"] * (th_des - thm) - G["kdt"] * wm)

        # actuator noise
        U = U * (1.0 + s["a_U"] * rng.standard_normal(n))
        tau = tau + s["a_tau"] * rng.standard_normal(n)

        # wind
        dx = p["wind"][0] + s["wind_gust"] * rng.standard_normal(n)
        dz = p["wind"][1] + s["wind_gust"] * rng.standard_normal(n)

        # true plant (random parameters)
        ax = -(U / p["m"]) * np.sin(th) + dx / p["m"]
        az = (U / p["m"]) * np.cos(th) - p["g"] + dz / p["m"]
        aw = tau / p["J"]

        # semi-implicit Euler
        vx, vz, w = vx + ax * DT, vz + az * DT, w + aw * DT
        x, z, th = x + vx * DT, z + vz * DT, th + w * DT
        if record:
            hist[k] = (x, z, th)

    err = np.hypot(x - xr, z - zr)
    return err, hist


def ci_mean(e):
    n = e.size
    half = Z_95 * e.std(ddof=1) / math.sqrt(n)
    return e.mean(), half


def ci_prop(success):
    n = success.size
    p = success.mean()
    return p, Z_95 * math.sqrt(p * (1 - p) / n)


def main():
    rng = np.random.default_rng(2026)

    print("STEP 1 - PVTOL model")
    print(f"  nominal parameters: {NOMINAL}")
    print(f"  controller gains:   {GAINS}")
    print(f"  target={TARGET}, dt={DT}s, T={T_FINAL}s, success if error < {SUCCESS_TOL} m")

    print_variance_table()

    print("\nSTEP 3 - Nominal simulation (all sigmas = 0)")
    e_nom, hist_nom = simulate(1, rng, scale=0.0, record=True)
    print(f"  final position error = {e_nom[0]:.6f} m")

    n0 = 500
    print(f"\nSTEP 4 - Pilot Monte Carlo, n0 = {n0}")
    e0, _ = simulate(n0, rng)
    s0 = e0.std(ddof=1)
    p0 = (e0 < SUCCESS_TOL).mean()
    print(f"  mean error = {e0.mean():.5f} m, sigma = {s0:.5f} m, p_success = {p0:.4f}")

    print(f"\nSTEP 5 - Minimum sample size (95 % confidence, E = {E_MAX})")
    n_mean = math.ceil((Z_95 * s0 / E_MAX) ** 2)
    n_prop = math.ceil(Z_95 ** 2 * p0 * (1 - p0) / E_MAX ** 2)
    n_worst = math.ceil(Z_95 ** 2 * 0.25 / E_MAX ** 2)
    print(f"  mean error : n = (z*sigma/E)^2      = ({Z_95:.3f}*{s0:.4f}/{E_MAX})^2 = {n_mean}")
    print(f"  proportion : n = z^2*p(1-p)/E^2     = {Z_95:.3f}^2*{p0:.3f}*{1 - p0:.3f}/{E_MAX}^2 = {n_prop}")
    print(f"  proportion worst case (p = 0.5)     = {n_worst}")
    n_req = max(n_mean, n_prop, n0)
    print(f"  -> N = max(n_mean, n_prop) = {n_req}")

    print(f"\nSTEP 6 - Full Monte Carlo, N = {n_req}")
    e, _ = simulate(n_req, rng)
    # the pilot sigma/p are estimates; top up if the final half-width still exceeds E
    while True:
        m, hm = ci_mean(e)
        p, hp = ci_prop(e < SUCCESS_TOL)
        if hm < E_MAX and hp < E_MAX:
            break
        n_more = max(math.ceil((Z_95 * e.std(ddof=1) / E_MAX) ** 2),
                     math.ceil(Z_95 ** 2 * p * (1 - p) / E_MAX ** 2)) - e.size + 1
        print(f"  half-width above E, adding {n_more} runs")
        e = np.concatenate([e, simulate(n_more, rng)[0]])

    print("\nSTEP 7 - Results")
    print(f"  N used               = {e.size}")
    print(f"  mean error           = {m:.5f} m, 95% CI [{m - hm:.5f}, {m + hm:.5f}], half-width = {hm:.5f}")
    print(f"  success probability  = {p:.4f}, 95% CI [{p - hp:.4f}, {p + hp:.4f}], half-width = {hp:.5f}")
    print(f"  both half-widths < {E_MAX}: {hm < E_MAX and hp < E_MAX}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    t = np.arange(hist_nom.shape[0]) * DT
    _, hist_mc = simulate(30, np.random.default_rng(1), record=True)
    i_all = np.arange(1, e.size + 1)
    running = np.cumsum(e) / i_all
    var = (np.cumsum(e ** 2) - i_all * running ** 2)[1:] / (i_all[1:] - 1)
    idx = i_all[1:]
    band = Z_95 * np.sqrt(np.maximum(var, 0.0)) / np.sqrt(idx)

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].plot(hist_mc[:, 0, :], hist_mc[:, 1, :], color="gray", alpha=0.4, lw=0.8)
    ax[0, 0].plot(hist_nom[:, 0, 0], hist_nom[:, 1, 0], "b", lw=2, label="nominal")
    ax[0, 0].plot(*TARGET, "r*", ms=12, label="target")
    ax[0, 0].set(xlabel="x [m]", ylabel="z [m]", title="Trajectories (30 MC runs)")
    ax[0, 0].legend()
    ax[0, 1].plot(t, hist_mc[:, 2, :], color="gray", alpha=0.4, lw=0.8)
    ax[0, 1].plot(t, hist_nom[:, 2, 0], "b", lw=2)
    ax[0, 1].set(xlabel="t [s]", ylabel="theta [rad]", title="Pitch angle")
    ax[1, 0].hist(e, bins=60, color="steelblue")
    ax[1, 0].axvline(SUCCESS_TOL, color="r", ls="--", label="success tol")
    ax[1, 0].set(xlabel="final position error [m]", ylabel="count", title=f"Error histogram, N={e.size}")
    ax[1, 0].legend()
    ax[1, 1].plot(idx, running[1:], "b", label="running mean")
    ax[1, 1].fill_between(idx, running[1:] - band, running[1:] + band, alpha=0.3, label="95% CI")
    ax[1, 1].axvline(n_req, color="r", ls="--", label="N required")
    ax[1, 1].set(xscale="log", xlabel="samples", ylabel="mean error [m]", title="Convergence")
    ax[1, 1].legend()
    fig.tight_layout()
    out = __file__.replace(".py", ".png")
    fig.savefig(out, dpi=120)
    print(f"\n  figure saved to {out}")


if __name__ == "__main__":
    main()

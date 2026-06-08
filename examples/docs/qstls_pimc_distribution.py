import numpy as np
import matplotlib.pyplot as plt

A_CUTOFF = 0.5


def load_pimc_data(rs: int):
    filename = f"Momentum_UEG_rs{rs}_theta1.txt"
    # File format: x y ylow yhigh type ; keep reliable "i" points when present.
    rows = []
    rows_i = []
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            x = float(parts[0])
            y = max(float(parts[1]), 0.0)
            rows.append((x, y))
            if len(parts) >= 5 and parts[4].lower() == "i":
                rows_i.append((x, y))
    data = np.array(rows_i if len(rows_i) >= 4 else rows, dtype=float)
    return data[:, 0], data[:, 1]


def g0_on_top(rs: float, theta: float = 1.0) -> float:
    theta_32 = theta**1.5
    theta2 = theta * theta
    theta3 = theta * theta2
    sqrt_theta = np.sqrt(theta)
    sqrt_rs = np.sqrt(rs)
    rs2 = rs * rs
    rs3 = rs2 * rs
    g0aa1 = 18.4376509088802
    g0ba1 = 24.1338558554951
    g0ba2 = 1.86499223116244
    g0a0 = 0.18315
    g0ab1 = -0.243679875065179
    g0bb1 = 0.252577
    g0bb2 = 0.12704315679703
    g0b0 = -0.0784043
    g0ac1 = 2.23662699250646
    g0ac2 = 0.448936594134834
    g0bc1 = 0.445525554738623
    g0bc2 = 0.408503601408896
    g0c0 = 1.02232
    g0ad1 = 0.0589015402409623
    g0bd1 = -0.598507865444083
    g0bd2 = 0.513161640681146
    g0d0 = 0.0837741
    g0a = (g0a0 + g0aa1 * theta) / (1.0 + g0ba1 * theta + g0ba2 * theta3)
    g0b = (g0b0 + g0ab1 * sqrt_theta) / (1.0 + g0bb1 * theta + g0bb2 * theta2)
    g0c = (g0c0 + g0ac1 * sqrt_theta + g0ac2 * theta_32) / (
        1.0 + g0bc1 * theta + g0bc2 * theta2
    )
    g0d = (g0d0 + g0ad1 * sqrt_theta) / (1.0 + g0bd1 * theta + g0bd2 * theta2)
    return 0.5 * (1.0 + g0a * sqrt_rs + g0b * rs) / (1.0 + g0c * rs + g0d * rs3)


def build_f0(y_data: np.ndarray, f_data: np.ndarray, rs: float, eta: float):
    lam = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
    axis_scale = lam * rs
    y_data = axis_scale * y_data
    y_sec = y_data[-1]
    g0 = g0_on_top(rs)
    tail_coeff = 4.0 / (9.0 * np.pi * np.pi) * lam * lam * rs * rs * g0
    y_max = y_data[-1]

    def f_pimc(y):
        return np.interp(np.clip(y, y_data[0], y_data[-1]), y_data, f_data)

    def f_tail(y):
        y = np.maximum(y, 1e-8)
        return tail_coeff / y**8

    def moment(y_sec_local):
        y_upper = max(100.0, 8.0 * y_max)
        y = np.linspace(0.0, y_upper, 6000)
        m = np.trapz(y * y * f0(y), y)
        m += tail_coeff / (5.0 * y_upper**5)
        return m

    def A(y):
        y = np.asarray(y, dtype=float)
        a = 0.5 * (1.0 + np.tanh(eta * (y - y_sec)))
        return np.where(y <= A_CUTOFF, 0.0, a)

    def f0(y):
        ay = A(y)
        return (1.0 - ay) * f_pimc(y) + ay * f_tail(y)

    # Solve y_sec from normalization with fixed eta.
    target = 1.0 / 3.0
    a, b = 0.01 * y_max, 10.0 * y_max
    def f_norm(yc):
        nonlocal y_sec
        y_sec = yc
        return moment(y_sec) - target
    fa, fb = f_norm(a), f_norm(b)
    if fa * fb > 0:
        scan = np.linspace(a, b, 220)
        vals = np.array([f_norm(x) for x in scan])
        idx = np.where(vals[:-1] * vals[1:] <= 0.0)[0]
        if idx.size > 0:
            i = int(idx[0])
            a, b = float(scan[i]), float(scan[i + 1])
            fa, fb = float(vals[i]), float(vals[i + 1])
        else:
            y_sec = float(scan[int(np.argmin(np.abs(vals)))])
            return y_sec, A, f0, moment(y_sec)
    for _ in range(70):
        c = 0.5 * (a + b)
        fc = f_norm(c)
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    y_sec = 0.5 * (a + b)
    return y_sec, A, f0, moment(y_sec)


def fermi_dirac_half_integral(mu: float) -> float:
    # Match native 3D normalization: Gamma(3/2) F_{1/2}(mu)
    t_max = max(60.0, mu + 60.0)
    t = np.linspace(0.0, t_max, 12000)
    integrand = np.sqrt(t) / (np.exp(t - mu) + 1.0)
    return np.trapz(integrand, t)


def chemical_potential_3d(theta: float = 1.0) -> float:
    target = 2.0 / (3.0 * theta ** 1.5)
    a, b = -20.0, 20.0
    fa = fermi_dirac_half_integral(a) - target
    fb = fermi_dirac_half_integral(b) - target
    while fa * fb > 0.0:
        a -= 10.0
        b += 10.0
        fa = fermi_dirac_half_integral(a) - target
        fb = fermi_dirac_half_integral(b) - target
    for _ in range(80):
        c = 0.5 * (a + b)
        fc = fermi_dirac_half_integral(c) - target
        if fa * fc <= 0.0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


def fermi_dirac_distribution(y: np.ndarray, theta: float = 1.0) -> np.ndarray:
    mu = chemical_potential_3d(theta)
    return 1.0 / (np.exp(y * y / theta - mu) + 1.0)


def main():
    rs = 10
    y_data, f_data = load_pimc_data(rs)
    etas = [4.0, 8.0, 16.0]
    lam = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
    axis_scale = lam * rs
    y_plot = np.linspace(0.0, max(8.0, 4.0 * axis_scale * y_data[-1]), 1200)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for eta in etas:
        y_sec, A, f0, norm = build_f0(y_data, f_data, rs, eta)
        print(f"eta={eta:.2f}, y_sec={y_sec:.4f}, moment={norm:.8f}")
        ax1.plot(
            y_plot,
            A(y_plot),
            label=fr"$\eta={eta:.1f}, y_{{sec}}={y_sec:.2f}$",
        )
        ax2.plot(
            y_plot,
            f0(y_plot),
            lw=2.0,
        label=fr"produced $f_0(y)$, $\eta={eta:.1f}$",
        )

    ax2.scatter(axis_scale * y_data, f_data, s=18, alpha=0.65, label="PIMC data")
    ax2.plot(
        y_plot,
        fermi_dirac_distribution(y_plot, 1.0),
        "--",
        lw=1.8,
        label="Fermi-Dirac",
    )

    ax1.set_xlabel("y")
    ax1.set_ylabel("A(y)")
    ax1.set_title("Switching function")
    ax1.grid(alpha=0.3)
    ax1.legend()

    ax2.set_xlabel("y")
    ax2.set_ylabel(r"$f_0(y)$")
    ax2.set_title(fr"Interacting distribution at $r_s={rs}, \theta=1$")
    ax2.grid(alpha=0.3)
    ax2.set_yscale("log")
    ax2.legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

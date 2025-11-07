import numpy as np
import matplotlib.pyplot as plt

# Constants
n_steps = 50000
step_size = 1.0
tol = 1e-4
a = 0.1
b = 1.0
lambda_ = 8.0  


def r12(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


def psi_2d(p1, p2, alpha):
    r = r12(p1, p2)
    gauss = np.exp(-0.5 * (np.dot(p1, p1) + np.dot(p2, p2)))
    jastrow = np.exp(lambda_ * r / (1 + alpha * r))
    return gauss * jastrow


def local_energy_2d(p1, p2, alpha):
    x1, y1 = p1
    x2, y2 = p2
    r = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    if r == 0:
        return np.inf 

    R2 = x1**2 + y1**2 + x2**2 + y2**2
    term1 = -R2
    term2 = -4
    term3 = 4*lambda_ / (((1 + alpha * r)**2)*r)
    term4 = -2 * ( ( ( lambda_ * (x1 - x2)**2) / (r**3 * (1 + alpha * r)**3) )* (1 + 3*alpha*r) )
    term5 = -2 * ( ( ( lambda_ * (y1 - y2)**2) / (r**3 * (1 + alpha * r)**3) )* (1 + 3*alpha*r) )
    term6 = ( -x1 + (lambda_ * (x1 - x2)) / ( (1 + alpha*r)**2 * r) )**2
    term7 = ( -x2 + (lambda_ * (x2 - x1)) / ( (1 + alpha*r)**2 * r) )**2
    term8 = ( -y1 + (lambda_ * (y1 - y2)) / ( (1 + alpha*r)**2 * r) )**2
    term9 = ( -y2 + (lambda_ * (y2 - y1)) / ( (1 + alpha*r)**2 * r) )**2
    
    E_tot = (-1/2) * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9) + lambda_ / r
    return E_tot



def metropolis_2d(alpha, n_steps, step_size):
    p1 = np.random.uniform(-1, 1, 2)
    p2 = np.random.uniform(-1, 1, 2)

    energies = []

    for _ in range(n_steps):
        new_p1 = p1 + step_size * np.random.uniform(-1, 1, 2)
        new_p2 = p2 + step_size * np.random.uniform(-1, 1, 2)

        p_old = psi_2d(p1, p2, alpha)**2
        p_new = psi_2d(new_p1, new_p2, alpha)**2
        acceptance = p_new / p_old

        if acceptance >= 1 or np.random.rand() < acceptance:
            p1, p2 = new_p1, new_p2

        energies.append(local_energy_2d(p1, p2, alpha))

    return np.array(energies)


def average_energy_2d(alpha):
    energies = metropolis_2d(alpha, n_steps, step_size)
    return np.mean(energies)


def golden_search(f, a, b, tol):
    gr = (np.sqrt(5) + 1) / 2
    c = b - (b - a) / gr
    d = a + (b - a) / gr

    while abs(b - a) > tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (b + a) / 2


def plot_energy_and_variance_2d():
    alphas = np.linspace(a, b, 100)
    optimal_alpha = golden_search(average_energy_2d, a, b, tol)
    print(f"Optimal alpha: {optimal_alpha:.4f}")

    energies = []
    variances = []

    for alpha in alphas:
        e = metropolis_2d(alpha, n_steps, step_size)
        energies.append(np.mean(e))
        variances.append(np.var(e))

    E_alpha = np.min(energies)
    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.set_xlabel("α")
    ax1.set_ylabel("Local energy", color="blue")
    ax1.plot(alphas, energies, color="blue", marker='o', label="⟨E_L⟩")
    ax1.axhline(E_alpha, color="gray", linestyle="--", label=f"Exact E = {E_alpha}")
    ax1.tick_params(axis='y', labelcolor="blue")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Variance", color="orange")
    ax2.plot(alphas, variances, color="orange", linestyle='--', marker='s', label="Variance")
    ax2.tick_params(axis='y', labelcolor="orange")

    ax1.axvline(optimal_alpha, color="red", linestyle=":", label=f"Optimal α ≈ {optimal_alpha:.3f}")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper center")

    plt.title(f"Energy and Variance for the 2D Harmonic oscillator (λ = {lambda_})")
    plt.tight_layout()
    plt.show()


def main():
    plot_energy_and_variance_2d()

if __name__ == "__main__":
    main()
    
    


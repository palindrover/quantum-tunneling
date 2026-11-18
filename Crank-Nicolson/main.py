import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

hbar = 1.0
m = 1.0
a = 1.0

def countPsi(V0, x):
    omega = np.sqrt(8 * V0 / (m * a**2))
    sigma = np.sqrt(hbar / (m * omega))

    dx = x[1] - x[0]

    psi = np.exp(-(x+a)**2 / (2 * sigma**2)).astype(complex)
    psi[0] = psi[-1] = 0

    C = 1 / np.sqrt(dx * np.sum(np.abs(psi)**2))
    psi *= C

    return psi


def main(V0, L, Nx, Tfin, Nt, n_out):
    save_every = max(1, Nt//n_out)

    x = np.linspace(-L, L, Nx)
    dx = x[1] - x[0]
    dt = Tfin/Nt

    V = V0 * (1 - x**2/a**2)**2

#    plt.plot(x, V)
#    plt.xlabel("x")
#    plt.ylabel("V(x)")
#    plt.show()

    r = 1j * hbar * dt / (4 * m * dx**2)
    beta = 1j * dt / (2 * hbar)

    Amain = 1 + 2*r + beta*V[1:-1]
    Asub = np.full(Nx-3, -r, dtype=complex)
    Bmain = 1 - 2*r - beta*V[1:-1]
    
    psi = countPsi(V0, x)

    ab = np.zeros((3, Nx - 2), dtype=complex)
    ab[0, 1:] = Asub
    ab[1] = Amain
    ab[2, :-1] = Asub

    time = []
    PL = []
    norm = []
    left = x < 0

    for i in range(Nt):
        if i % save_every == 0:
            time.append(i*dt)
            PL.append(dx * np.sum(np.abs(psi[left])**2))
            norm.append(dx * np.sum(np.abs(psi)**2))
        
        Bpsi = Bmain * psi[1:-1] + r * (psi[:-2] + psi[2:])
        psi[1:-1] = solve_banded((1, 1), ab, Bpsi)

    plt.plot(time, PL, label="$P_L(t)$")
    plt.plot(time, norm, label="norm")
    plt.legend()
    plt.grid(True)
    plt.show()

main(8.0, 3.0, 1001, 100.0, 10000, 2000)

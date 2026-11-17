import numpy as np
import matplotlib.pyplot as plt

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

    plt.plot(x, np.abs(psi)**2)
    plt.show()

    return psi


def main(V0, L, Nx, Tfin, Nt):
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

    Amain = 1 + 2*r + beta*V
    Asub = np.full(Nx-3, -r, dtype=complex)
    Bmain = 1 - 2*r - beta*V
    Bsub = np.full(Nx-3, r, dtype=complex)
    
    psi = countPsi(V0, x)

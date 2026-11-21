import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded, eigh_tridiagonal

hbar = 1.0
m = 1.0
a = 1.0

def graph(x, V, psi):
    plt.plot(x, V, label="V(x)")
    plt.plot(x, psi, label="$\psi(x)$")
    plt.plot(x, np.abs(psi)**2, label="$|\psi|^2$")
    plt.legend()
    plt.grid(True)
    plt.show()

def countPsi(V0, x):
    omega = np.sqrt(8 * V0 / (m * a**2))
    sigma = np.sqrt(hbar / (m * omega))

    dx = x[1] - x[0]

    psi = np.exp(-(x+a)**2 / (2 * sigma**2)).astype(complex)
    psi[0] = psi[-1] = 0

    C = 1 / np.sqrt(dx * np.sum(np.abs(psi)**2))
    psi *= C

    return psi

def eigenPsi(V, x):
    dx = x[1] - x[0]
    Hmain = hbar**2 / (m * dx**2) + V[1:-1]
    Hsub = np.full(len(x)-3, -1 * hbar**2 / (2 * m * dx**2))

    w, v = eigh_tridiagonal(Hmain, Hsub, select='i', select_range=(0, 1))

    psi = np.zeros(len(x), dtype=complex)
    psi[1:-1] = (v[:, 0] + v[:, 1]) / np.sqrt(2) 

    C = 1 / np.sqrt(dx * np.sum(np.abs(psi)**2))
    psi *= C

    left = x < 0
    if dx * np.sum(np.abs(psi[left])**2) < 0.5:
        psi = psi[::-1]

    return psi, w


def measurePeriod(time, PL):
    time = np.array(time)
    PL = np.array(PL)
    f = PL - 0.5
    idx = np.where(np.diff(np.sign(f)) != 0)[0]
    t = time[idx] - f[idx] * (time[idx+1] - time[idx]) / (f[idx+1] - f[idx]) 
    
    return np.mean(t[2:] - t[:-2])

def buildSystem(V0, L, Nx, Tfin, Nt):
    x = np.linspace(-L, L, Nx)
    dx = x[1] - x[0]
    dt = Tfin/Nt

    V = V0 * (1 - x**2/a**2)**2

    r = 1j * hbar * dt / (4 * m * dx**2)
    beta = 1j * dt / (2 * hbar)

    Amain = 1 + 2*r + beta*V[1:-1]
    Asub = np.full(Nx-3, -r, dtype=complex)
    Bmain = 1 - 2*r - beta*V[1:-1]

    ab = np.zeros((3, Nx - 2), dtype=complex)
    ab[0, 1:] = Asub
    ab[1] = Amain
    ab[2, :-1] = Asub

    return x, dx, dt, V, r, ab, Bmain

def evolve(psi, x, dx, dt, Nt, n_out, Bmain, ab, r):
    psi = psi.copy()
    save_every = max(1, Nt//n_out)
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

    return time, PL, norm

def compareSplitting(w, dt, T):
    dE = w[1] - w[0]
    dECN = hbar * 2 * (np.arctan(w[1] * dt / (2 * hbar)) - 
                       np.arctan(w[0] * dt / (2 * hbar))) / dt
    dEm = 2 * np.pi * hbar / T

    methodError = (dECN - dE) / dE
    measureError = (dEm - dECN) / dECN

    return {
            "dE": dE, 
            "dECN": dECN, 
            "dEm": dEm, 
            "methodError": methodError, 
            "measureError": measureError
            }

def main(V0, L, Nx, Tfin, Nt, n_out):

    x, dx, dt, V, r, ab, Bmain = buildSystem(V0, L, Nx, Tfin, Nt)
    
    psi, w = eigenPsi(V, x)

    #if Nt == 20000: graph(x, V, psi)

    time, PL, norm = evolve(psi, x, dx, dt, Nt, n_out, Bmain, ab, r)

    T = measurePeriod(time, PL)

    res = compareSplitting(w, dt, T)

    print(f"\n\n{Nt}\n")
    for key, val in res.items():
        print(f"{key:12s} {val: .6e}")

if __name__ == "__main__":
    for Nt in [2000, 5000, 10000, 20000]:
        main(8.0, 3.0, 1001, 100.0, Nt, 2000)

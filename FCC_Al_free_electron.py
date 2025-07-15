import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Publication-style settings
rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.linewidth': 1.2,
    'lines.linewidth': 1.2,
    'xtick.major.width': 1.1,
    'ytick.major.width': 1.1,
    'figure.figsize': (6, 5),
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'text.usetex': False
})

# Physical constants
hbar = 1 
m = 1   

# FCC lattice constant Aluminum
a_angstrom = 4.05  # Å
a = a_angstrom # meters

vol = a_angstrom**3/4
scaling = vol**(float(-2/3))

# Reciprocal lattice vectors in 1/m
b1 = 2 * np.pi / a * np.array([1, 1, -1])
b2 = 2 * np.pi / a * np.array([1, -1, 1])
b3 = 2 * np.pi / a * np.array([-1, 1, 1])

# High-symmetry points
G = np.array([0, 0, 0])
K = np.array([3*np.pi/(2*a),3*np.pi/(2*a),0])

# Interpolate k-path from Γ to K
n_kpoints = 200
k_path = np.linspace(G, K, n_kpoints)

# Reciprocal lattice vectors for band folding
G_vectors = []
max_index = 1
for i in range(-max_index, max_index + 1):
    for j in range(-max_index, max_index + 1):
        for k in range(-max_index, max_index + 1):
            G_vec = i*b1 + j*b2 + k*b3
            G_vectors.append(G_vec)

# Calculate energies in eV
E_bands = []
for G_vec in G_vectors:
    k_plus_G = k_path + G_vec
    k2 = np.sum(k_plus_G**2, axis=1)
    E = (hbar**2 / (2 * m)) * k2
    E_bands.append(E)

# Plot
fig, ax = plt.subplots()
for band in E_bands:
    ax.plot(np.linspace(0, 1, n_kpoints), band*scaling, color='navy', lw=1)

ax.set_title(r'Free Electron Bands (FCC)')
ax.set_ylabel(r'Energy (Rydberg $\cdot \Omega^{-2/3}$)')
plt.ylim(0,1)
ax.set_xticks([0, 1])
ax.set_xticklabels([r'$\Gamma$', r'$K$'])
ax.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig("gamma_to_k.png")
plt.show()

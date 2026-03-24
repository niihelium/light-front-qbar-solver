import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib.pyplot as plt

def plot_eigenstates(evals, evecs, r_axis, w_axis, resolution, num_states):
    """
    Plots the first `num_states` eigenstates in a grid layout.
    """
    # Grid layout: 4 columns, rows calculated based on num_states
    cols = 4
    rows = (num_states + cols - 1) // cols
    
    # Increased figsize for better resolution and readability
    fig, axes = plt.subplots(rows, cols, figsize=(20, 5 * rows), sharex=True, sharey=True, squeeze=False)
    
    # Adjust subplots to give space for larger labels
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.12, wspace=0.15, hspace=0.3)

    for i in range(num_states):
        ax = axes[i // cols, i % cols]
        
        # Extract 1D vector and reshape into 2D grid (indices: [r, w])
        # Transpose (.T) so that r is horizontal (X) and w is vertical (Y)
        state_2d = evecs[:, i].reshape((resolution, resolution)).T
        
        # Plotting
        ax.imshow(state_2d, origin='lower', extent=[r_axis[0], r_axis[-1], w_axis[0], w_axis[-1]], aspect='auto', cmap='viridis')
        
        # Set title with energy value
        ax.set_title(f"State {i+1}, E = {evals[i]:.2f}", fontsize=18, fontweight='bold', pad=10)
        
        # Increase tick label size
        ax.tick_params(axis='both', which='major', labelsize=14)

    # Label outer axes
    for ax in axes[-1, :]:
        ax.set_xlabel('r', fontsize=22, labelpad=10)
    for ax in axes[:, 0]:
        ax.set_ylabel('w', fontsize=22, labelpad=10)

    # Hide any unused subplots
    for j in range(num_states, rows * cols):
        axes[j // cols, j % cols].axis('off')

    # Save as image with tight layout
    filename = 'energy_states.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Results saved to {filename}")

# === 1. Main Parameters ===
num_states = 8
resolution = 128

# Boundaries 
w_interval = (-1.0, 1.0)
r_interval = (0.0, 1.0)

# === 2. Math Calculations ===

# Axis w (using standard linspace)
w_axis, w_step = np.linspace(w_interval[0], w_interval[1], resolution, retstep=True)

# Axis r (offset from r=0 by one grid step since u(0) = 0)
r_step = (r_interval[1] - r_interval[0]) / (resolution + 1)
r_axis = np.linspace(r_interval[0] + r_step, r_interval[1] - r_step, resolution)

# Constructing 1D Operators
# Radial axis (r)
D2_r = sp.diags([1, -2, 1], [-1, 0, 1], shape=(resolution, resolution)) / (r_step**2)
inv_r2 = sp.diags(1.0 / (r_axis**2))
V_r = sp.diags(r_axis**2)  # Harmonic potential along r

# Symmetrized Hamiltonian
# Substitution: psi = u / sqrt(r) -> operator becomes -D^2_r - 1/(4*r^2)
H_r_sym = -D2_r - 0.25 * inv_r2 + V_r

# Vertical axis (w)
D2_w = sp.diags([1, -2, 1], [-1, 0, 1], shape=(resolution, resolution)) / (w_step**2)
V_w = sp.diags(w_axis**2)
H_w = -D2_w + V_w

# Assembling Full 2D Hamiltonian
# Kronecker product for 2D grid: H_2d = H_r (x) I_w + I_r (x) H_w
I_res = sp.eye(resolution)
H_2d = sp.kron(H_r_sym, I_res) + sp.kron(I_res, H_w)

# Finding Eigenvalues
print(f"Finding {num_states} states for Hermitian matrix of size {H_2d.shape}...")
evals, evecs = sla.eigsh(H_2d, k=num_states, which='SM')

# Sort levels by energy (ascending order)
idx = np.argsort(evals)
evals = evals[idx]
evecs = evecs[:, idx]

# === 3. Visualization ===
plot_eigenstates(evals, evecs, r_axis, w_axis, resolution, num_states)
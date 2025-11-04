import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# Constants from Fibonacci Cosmology Theorem
G = 1.0  # Gravitational constant (normalized units)
phi = (1 + np.sqrt(5)) / 2  # Golden Ratio
phi_conj = phi - 1  # Conjugate (1/phi)
t0 = 1.0  # Characteristic time (normalized; in real units ~1/H0)
rho_crit = 0.5  # Critical density threshold for phase flip (tunable)
dt = 0.01  # Time step
n_steps = 500  # Number of simulation steps
n_particles = 50  # Number of particles
box_size = 10.0  # Initial box size for rough density normalization
epsilon = 0.1  # Softening length to avoid force singularities

# Initialize particles: positions, velocities, masses
np.random.seed(42)
pos = np.random.uniform(0, box_size, (n_particles, 2))
vel = np.random.normal(0, 0.1, (n_particles, 2))
mass = np.ones(n_particles)  # Unit masses for simplicity


# Function to compute softened gravitational forces
def compute_forces(pos, mass):
    dist_vec = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
    dist_sq = np.sum(dist_vec ** 2, axis=2) + epsilon ** 2  # Softened
    dist_sq[np.diag_indices(n_particles)] = np.inf
    force_mag = G * (mass[:, np.newaxis] * mass) / dist_sq
    forces = np.sum(force_mag[:, :, np.newaxis] * (dist_vec / np.sqrt(dist_sq)[:, :, np.newaxis]), axis=1)
    return forces


# Function to estimate local density (average inverse distance squared proxy)
def estimate_density(pos):
    dist = pdist(pos)
    dist[dist == 0] = np.inf
    return np.mean(1 / dist ** 2)  # Improved proxy for density


# Simulation loop
positions_over_time = [pos.copy()]
a_over_time = [1.0]  # Initial scale factor
sigma_over_time = [1.0]  # Initial phase (expansion)

for step in range(n_steps):
    # Compute current density
    rho = estimate_density(pos)

    # Determine phase sigma (smooth transition for stability: tanh)
    sigma = np.tanh(rho_crit - rho)  # ~ +1 if rho < crit, -1 if rho > crit
    sigma_over_time.append(sigma)

    # Compute scale factor change based on phase
    ln_r = np.log(phi)  # Magnitude same for dual (approx)
    da_dt = sigma * (ln_r / t0)
    a_new = a_over_time[-1] + da_dt * dt
    a_over_time.append(a_new)

    # Compute forces
    forces = compute_forces(pos, mass)

    # Update velocities (Leapfrog integrator + Hubble term)
    H = da_dt / a_over_time[-1]  # Hubble parameter
    vel += (forces / mass[:, np.newaxis]) * dt / 2  # Half-kick
    pos += vel * dt
    forces = compute_forces(pos, mass)  # Recompute after drift
    vel += (forces / mass[:, np.newaxis]) * dt / 2 - H * vel * dt  # Full kick + drag

    # Apply cosmological scaling to positions
    scale_ratio = a_new / a_over_time[-2]
    pos *= scale_ratio

    positions_over_time.append(pos.copy())

# Plot results for visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Particle trajectories
for i in range(n_particles):
    traj = np.array([p[i] for p in positions_over_time])
    ax1.plot(traj[:, 0], traj[:, 1], alpha=0.5)
ax1.scatter(pos[:, 0], pos[:, 1], c='red', label='Final Positions')
ax1.set_title('Particle Trajectories in Dual-Phase Patch')
ax1.set_xlabel('x');
ax1.set_ylabel('y')
ax1.legend()
ax1.grid(True)

# Scale factor and phase evolution
ax2.plot(a_over_time, label='Scale Factor a(t)')
ax2.plot(sigma_over_time, '--', label='Phase σ(t)')
ax2.set_title('Evolution of Scale Factor and Phase')
ax2.set_xlabel('Time Step');
ax2.set_ylabel('Value')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('n_body_fib_dual.png', dpi=150)
print("Simulation complete. Plot saved as 'n_body_fib_dual.png'")

# Summary statistics
print(f"Final density: {estimate_density(pos):.3f}")
print(f"Effective phase flips (crossings): {sum(np.diff(np.sign(sigma_over_time)) != 0)}")

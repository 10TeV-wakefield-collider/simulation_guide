import numpy as np
import scipy.constants as sc
from wake_t import ParticleBunch, GaussianPulse, PlasmaStage, Beamline

# Plasma density profile
# ----------------------
n_p0 = 1.7e23 # On-axis density in the plasma channel (m^-3)
kp = np.sqrt(n_p0 * sc.e**2 / (sc.epsilon_0 * sc.m_e * sc.c**2))  # plasma wavenumber
l_plateau = 28e-2  # plateau length
l_ramp = 1e-3  # ramp length
l_total = l_plateau + 2 * l_ramp
# Determine guiding channel.
r_e = sc.e**2 / (4.0 * np.pi * sc.epsilon_0 * sc.m_e * sc.c**2)
matched_channel_waist = 40.e-6 # determined empirically for best laser guiding
# matched_channel_waist should be equal to w_0 for low a0, but differs for high a0 due to relativistic guiding
rel_delta_n_over_w2 = 1.0 / (np.pi * r_e * matched_channel_waist ** 4 * n_p0)

# Density function:
# A simple flat-top plasma density region with cosine ramps and a parabolic radial profile.
def density_profile(z, r):
    """ Define plasma density as a function of ``z`` and ``r``. """
    # Allocate density array.
    n = n_p0 * np.ones_like(z)
    # Up ramp
    n = np.where((z >= 0.0) & (z <= l_ramp), 
                 0.5 * (1 - np.cos(np.pi * z / l_ramp)) * n, n)
    # Down ramp
    z0 = l_ramp + l_plateau
    n = np.where((z >= z0) & (z <= z0 + l_ramp), 
                 (1 - 0.5 * (1 - np.cos(np.pi * (z - z0) / l_ramp))) * n, n)
    # Set to a very small, non-zero value outside the plasma border
    n = np.where(z <= 0.0, 1e-10 * n_p0, n)
    n = np.where(z >= z0 + l_ramp, 1e-10 * n_p0, n)
    n = np.where(n == 0, 1e-10 * n_p0, n)  # fix having zero in the last point
    # Add radial parabolic profile
    n = n * (1. + rel_delta_n_over_w2 * r**2)
    return n

# Laser setup
# -----------
a_0 = 2.36
w_0 = 36e-6
tau_fwhm = 7.33841e-14 * np.sqrt(2.0 * np.log(2))  # FWHM in intensity
xi_c = 0  # Center of the drive laser pulse
lambda_L = 0.8e-6  # laser wavelength (m)
laser = GaussianPulse(xi_c=xi_c, a_0=a_0, w_0=w_0,
                      l_0=lambda_L, tau=tau_fwhm,
                      z_foc=l_ramp)

# Witness beam setup
# ------------------
q_tot = 196e-12  # total charge (not counting the gaussian flanks)
ene_sp = 0.1  # energy spread in %
n_emitt_x = 0.1e-6
n_emitt_y = 0.1e-6

kin_energy_GeV = 100.0  # initial beam kinetic energy
mc2_GeV = sc.m_e * sc.c**2 / sc.electron_volt / 1e9
gamma_beam = 1 + kin_energy_GeV / mc2_GeV
betax0 = np.sqrt(2 * gamma_beam) / kp  # matched beta for a perfect blowout
sx0 = np.sqrt(n_emitt_x * betax0 / gamma_beam)  # matched beam size (rms)
sy0 = np.sqrt(n_emitt_y * betax0 / gamma_beam)  # matched beam size (rms)

# Best values from the case of 086_APOSSMM..
l_beam = 5.453219e-6  # beam length (m)
sz_beam = 1e-7  # gaussian flank sigma
d_beam = 62.973740e-6  # distance from the drive laser center to the beam front
i_r1 = 0.53034  # frontal current (relative)
i_r0 = 1.0 - i_r1  # rear current (relative)
i_avg = q_tot * sc.c / l_beam  # average beam current
i1_beam = 2 * i_r1 * i_avg  # frontal current (A)
i0_beam = 2 * i_r0 * i_avg  # rear current (A)

n_part = 50000

# Fix random seed for reproducibility
np.random.seed(0)

def trapezoidal_bunch( i0, i1, n_part, gamma0, s_g, length, s_z, emit_x, s_x, emit_y, s_y):
    """Create a trapezoidal particle bunch."""
    n_part = int(n_part)

    q_plat = (min(i0, i1) / sc.c) * length
    q_triag = ((max(i0, i1) - min(i0, i1)) / sc.c) * length / 2.0
    q_gaus0 = (i0 / sc.c) * np.sqrt(2 * np.pi) * s_z / 2.0
    q_gaus1 = (i1 / sc.c) * np.sqrt(2 * np.pi) * s_z / 2.0
    q_tot = q_plat + q_triag + q_gaus0 + q_gaus1

    n_plat = int(n_part * q_plat / q_tot)
    n_triag = int(n_part * q_triag / q_tot)
    n_gaus0 = int(n_part * q_gaus0 / q_tot)
    n_gaus1 = int(n_part * q_gaus1 / q_tot)

    z_plat = np.random.uniform(0.0, length, n_plat)
    if i0 <= i1:
        z_triag = np.random.triangular(0., length, length, n_triag)
    else:
        z_triag = np.random.triangular(0., 0., length, n_triag)
    z_gaus0 = s_z * np.random.standard_normal(2 * n_gaus0)
    z_gaus1 = s_z * np.random.standard_normal(2 * n_gaus1)

    z = np.concatenate((
        z_gaus0[np.where(z_gaus0 < 0)],
        z_plat, z_triag,
        z_gaus1[np.where(z_gaus1 > 0)] + length,
    ))
    z = z - length / 2.0  # shift to final position

    n_part = len(z)
    x = s_x * np.random.standard_normal(n_part)
    y = s_y * np.random.standard_normal(n_part)

    gamma = np.random.normal(gamma0, s_g, n_part)
    s_ux = emit_x / s_x
    ux = s_ux * np.random.standard_normal(n_part)
    s_uy = emit_y / s_y
    uy = s_uy * np.random.standard_normal(n_part)
    uz = np.sqrt((gamma**2 - 1) - ux**2 - uy**2)
    q = np.ones(n_part) * q_tot / n_part

    return x, y, z, ux, uy, uz, q

# Generate bunch particles
x, y, z, ux, uy, uz, q = trapezoidal_bunch(
    i0_beam,
    i1_beam,
    n_part=n_part,
    gamma0=gamma_beam,
    s_g=ene_sp * gamma_beam / 100,
    length=l_beam,
    s_z=sz_beam,
    emit_x=n_emitt_x,
    s_x=sx0,
    emit_y=n_emitt_y,
    s_y=sy0,
)
zc = xi_c - d_beam - l_beam / 2  # center position of the beam
zf = l_ramp  # beam focus position 
ctf = zf - zc  # distance to the focal plane
# backpropagate the beam balistically from the focal plane to the start of the simulation
gamma = np.sqrt(1 + ux**2 + uy**2 + uz**2)
z = z + zf - uz * ctf / gamma
x = x - ux * ctf / gamma
y = y - uy * ctf / gamma
w = np.abs(q / sc.e)
bunch = ParticleBunch(w, x, y, z, ux, uy, uz, name="bunch")

# Create plasma stage
# -------------------
# Simulation box and grid parameters.
r_max = w_0 * 4  # Maximum radial extent of the box
r_max_plasma = w_0 * 3  # Maximum radial extent of the plasma.
l_box = 200e-6  # Length of simulation box
xi_max = xi_c + 45e-6  # Right edge of the box in the speed-of-light frame
xi_min = xi_max - l_box  # Left edge of the box in the speed-of-light frame
dz = 1 / kp / 40
nz = int(l_box / dz)
dr = 1 / kp / 20
nr = int(r_max / dr)
dz_fields = 200e-6  # wakefield calculation period

# Adaptive grid setup
res_beam_r = 5.  # Beam radial resolution: number of grid points per sigma
adaptive_dr = np.min([sx0, sy0]) / res_beam_r
sxi = np.std(x)  # initial beam rms size in x
syi = np.std(y)  # initial beam rms size in y
adaptive_grid_r_max = 4 * np.max([sxi, syi])
adaptive_grid_nr = int(adaptive_grid_r_max / adaptive_dr)
# add more particles in the adaptive grid region to match the increased resolution
# Extends the region with more particles a bit beyond the adaptive grid border
ppc = 2
ppc = [[1.5 * adaptive_grid_r_max, ppc * int(dr / adaptive_dr)],
       [r_max_plasma, ppc]]

# Create plasma stage
plasma_plateau = PlasmaStage(
    length=l_total,
    density=density_profile,
    wakefield_model="quasistatic_2d",
    ion_motion=True,
    n_out=50,
    laser=laser,
    laser_evolution=True,
    r_max=r_max,
    r_max_plasma=r_max_plasma,
    xi_min=xi_min,
    xi_max=xi_max,
    n_r=nr,
    n_xi=nz,
    dz_fields= dz_fields,
    ppc=ppc,
    laser_envelope_substeps=4,
    laser_envelope_nxi=nz * 4,
    max_gamma=25,
    field_diags=["rho", "E", "B", "a"],
    use_adaptive_grids=True,
    adaptive_grid_r_max=adaptive_grid_r_max,
    adaptive_grid_nr=adaptive_grid_nr,
    adaptive_grid_diags=["E", "B"],
)

if __name__ == "__main__":

    # Track the beam through the plasma stage
    beamline = Beamline([ plasma_plateau ])
    beamline.track(bunch, opmd_diag=True, show_progress_bar=True, diag_dir='diags')

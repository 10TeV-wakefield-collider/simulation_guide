import numpy as np
import scipy.constants as sc
from wake_t import ParticleBunch, GaussianPulse, PlasmaStage, Beamline

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
    z_triag = np.random.triangular(0.0, length, length, n_triag)
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

def beta_gamma_from_GeV(m_species,E_kin_GeV):
    E_rest_GeV = m_species * sc.c**2 / sc.electron_volt / 1e9
    return np.sqrt((E_kin_GeV/E_rest_GeV + 1)**2 - 1)

def gamma_from_GeV(m_species,E_kin_GeV):
    E_rest_GeV = m_species * sc.c**2 / sc.electron_volt / 1e9
    return 1 + E_kin_GeV / E_rest_GeV

def calculate_critical_density(wavelength):
    """Calculate the critical density for a given laser wavelength."""
    omega_L = 2 * np.pi * sc.c / wavelength
    return sc.m_e * sc.epsilon_0 * omega_L**2 / sc.elementary_charge**2

def matched_beam_size(m_species, E_kin_GeV, beta_matched, norm_emittance):
    """ Matched beam size (1-RMS) for full blowout. """
    bg = beta_gamma_from_GeV(m_species, E_kin_GeV)
    return np.sqrt(norm_emittance * beta_matched / bg)

# General simulation parameters.
n_p0 = 1.7e23
q_tot = 196 * 1e-12  # total charge (C)
lambda_L = 0.8e-6  # laser wavelength (m)

# Energy values and corresponding beta function values from calculations and measurements
#   (previous optimization at 1 fC charge)
energies_GeV = np.array([1,10,100,1e3,1e4])  # initial beam energy
betas_fb = np.array([8.065e-4, 2.550e-3, 8.063e-3, 2.550e-2, 8.063e-2])  # analytical value for full blowout
betas_opt = np.array([9.330e-4, 2.775e-3, 7.746e-6, 1.721e-2, 3.644e-2])  # measured in optimization study
matched_beam_size(m_species=sc.m_e, E_kin_GeV=1, beta_matched=betas_opt[0], norm_emittance=1e-6)

# Best values from the case of 086_APOSSMM..
beam_i_r2 = 0.469660
beam_z_i_2 = 2.973740
beam_length = 5.453219

# Laser parameters
# ----------------
a_0 = 2.36
w_0 = 36e-6
tau_fwhm = 7.33841e-14 * np.sqrt(2.0 * np.log(2))  # FWHM in intensity
z_foc = 0.0  # Focal position of the laser

# Beam parameters
# ---------------
i = 0  # chooses the beam energy (1 GeV for i = 0) and beta values from above
z_beam = (60 + beam_z_i_2) * 1e-6  # distance of the beam center to the drive laser pulse
l_beam = beam_length * 1e-6  # beam length (m)

# beam currents are coupled due to total charge
beam_i_r1 = 1.0 - beam_i_r2

# we keep the beam charge almost constant (just not accounting for the gaussian flanks yet)
i_avg = q_tot / l_beam * sc.c  # average beam current (for rectangle profile)
i1_beam = 2 * beam_i_r1 * i_avg
i2_beam = 2 * beam_i_r2 * i_avg

n_part = 50000
gamma_beam = gamma_from_GeV(m_species=sc.electron_mass, E_kin_GeV=energies_GeV[i])  # calculate beam gamma
ene_sp = 0.1  # energy spread in %
sz_beam = 1e-7
n_emitt_x = 1e-6
n_emitt_y = 1e-6
betax0 = betas_opt[i]

sx0 = np.sqrt(n_emitt_x * betax0 / gamma_beam)  # matched beam size (rms)
sy0 = np.sqrt(n_emitt_y * betax0 / gamma_beam)  # matched beam size (rms)

# Generate bunch
x, y, z, ux, uy, uz, q = trapezoidal_bunch(
    i1_beam,
    i2_beam,
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

z -= l_beam / 2 + z_beam  # shift the longitudinal beam position away from the driver
w = np.abs(q / sc.e)
bunch = ParticleBunch(w, x, y, z, ux, uy, uz, name="bunch")

# Plasma density profile
# ----------------------
l_plateau = 28.0e-2 # (m)
# Determine guiding channel.
r_e = sc.e**2 / (4.0 * np.pi * sc.epsilon_0 * sc.m_e * sc.c**2)

# empirical length scale for guiding the pulse in this partial blowout
emp_len_um = 40.e-6
rel_delta_n_over_w2 = 1.0 / (np.pi * r_e * emp_len_um ** 4 * n_p0)

# Density function.
def density_profile(z):
    """ Define plasma density as a function of ``z``. """
    # Allocate relative density array.
    n = np.zeros_like(z)
    # Add plateau.
    n = np.where(
        (z >= 0) & (z <= l_plateau),
        1, n
    )
    # Return absolute density.
    n_abs = n * n_p0
    # Add minimal density value everywhere where it is zero, so we avoid crashes
    # Date: 2024-01-10, recommended by Angel Ferran-Pousa
    n_fixed = np.where(n_abs <= 1e10, 1e10, n_abs)
    return n_fixed

# Create plasma stage
# -------------------
# Simulation box and grid parameters.
r_max = w_0 * 4  # Maximum radial extent of the box
r_max_plasma = w_0 * 3  # Maximum radial extent of the plasma.
l_box = 200e-6  # Length of simulation box
# make sure initialize the bunch with a distance to the driver, or change this number
xi_0 = 0  # Center of the drive laser pulse
xi_max = xi_0 + 45e-6  # Right edge of the box in the speed-of-light frame
xi_min = xi_max - l_box  # Left edge of the box in the speed-of-light frame
res_beam_r = 5.  # default: 5, resolution we want for the beam radius
dr = np.min([sx0,sy0]) / res_beam_r  # make the radial resolution depend on the RMS beam size to avoid artifacts

n_e_n_c = n_p0 / calculate_critical_density(wavelength=lambda_L)
kp = 2 * np.pi * np.sqrt(n_e_n_c) / lambda_L
kpdz_inv = 40 # 1 / (kp * dz), resolution parameter where typically 20 - 40 is good
dz = 1 / kpdz_inv / kp
nr = int(r_max / dr)
nz = int(l_box / dz)
laser_evolution = True

laser = GaussianPulse(xi_c=xi_0, a_0=a_0, w_0=w_0, l_0=lambda_L, tau=tau_fwhm, z_foc=z_foc)

# Create plasma stages.
plasma_plateau = PlasmaStage(
    length=l_plateau,
    density=density_profile,
    wakefield_model="quasistatic_2d",
    n_out=50,
    laser=laser,
    laser_evolution=laser_evolution,
    r_max=r_max,
    r_max_plasma=r_max_plasma,
    xi_min=xi_min,
    xi_max=xi_max,
    n_r=nr,
    n_xi=nz,
    dz_fields= 1 * l_box,
    ppc=2,
    parabolic_coefficient=rel_delta_n_over_w2,
    laser_envelope_substeps=4,
    laser_envelope_nxi=nz * 4,
    max_gamma=25
)

# Track the beam through the plasma stage
beamline = Beamline([ plasma_plateau ])
beamline.track(bunch, opmd_diag=True, show_progress_bar=True, diag_dir='diags')

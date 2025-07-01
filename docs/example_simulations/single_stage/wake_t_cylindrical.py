import numpy as np
import scipy.constants as sc
from wake_t import ParticleBunch, GaussianPulse, PlasmaStage, Beamline

# Laser setup
# -----------
a_0 = 2.36
w_0 = 36e-6
tau_fwhm = 7.33841e-14 * np.sqrt(2.0 * np.log(2))  # FWHM in intensity
xi_0 = 0  # Center of the drive laser pulse
lambda_L = 0.8e-6  # laser wavelength (m)
laser = GaussianPulse(xi_c=xi_0, a_0=a_0, w_0=w_0, l_0=lambda_L, tau=tau_fwhm)

# Plasma density profile
# ----------------------
n_p0 = 1.7e23 # On-axis density in the plasma channel (m^-3)
kp = np.sqrt(n_p0 * sc.e**2 / (sc.epsilon_0 * sc.m_e * sc.c**2))  # plasma wavenumber
l_plateau = 28.e-2 # (m)
# Determine guiding channel.
r_e = sc.e**2 / (4.0 * np.pi * sc.epsilon_0 * sc.m_e * sc.c**2)
matched_channel_waist = 40.e-6 # determined empirically for best laser guiding
# matched_channel_waist should be equal to w_0 for low a0, but differs for high a0 due to relativistic guiding
rel_delta_n_over_w2 = 1.0 / (np.pi * r_e * matched_channel_waist ** 4 * n_p0)

# Density function:
# A simple constant plasma density with a parabolic radial profile.
def density_profile(z, r):
    """ Define plasma density as a function of ``z`` and ``r``. """
    # Allocate density array.
    n = n_p0 * np.ones_like(z)
    # Add radial parabolic profile
    n = n * (1. + rel_delta_n_over_w2 * r**2)
    return n

# Witness beam setup
# ------------------
q_tot = 196 * 1e-12  # total charge (C)
ene_sp = 0.1  # energy spread in %
n_emitt_x = 1e-6
n_emitt_y = 1e-6

kin_energy_GeV = 1.0  # initial beam kinetic energy
mc2_GeV = sc.m_e * sc.c**2 / sc.electron_volt / 1e9
gamma_beam = 1 + kin_energy_GeV / mc2_GeV
betax0 = np.sqrt(2 * gamma_beam) / kp  # matched beta for a perfect blowout
sx0 = np.sqrt(n_emitt_x * betax0 / gamma_beam)  # matched beam size (rms)
sy0 = np.sqrt(n_emitt_y * betax0 / gamma_beam)  # matched beam size (rms)

# Best values from the case of 086_APOSSMM..
l_beam = 5.453219e-6  # beam length (m)
beam_z_i_2 = 2.973740
z_beam = (60 + beam_z_i_2) * 1e-6  # distance of the beam center to the drive laser pulse
# beam currents are coupled due to total charge
beam_i_r2 = 0.469660
beam_i_r1 = 1.0 - beam_i_r2
sz_beam = 1e-7  # gaussian flank size
# we keep the beam charge almost constant (just not accounting for the gaussian flanks yet)
i_avg = q_tot / l_beam * sc.c  # average beam current (for rectangle profile)
i1_beam = 2 * beam_i_r1 * i_avg
i2_beam = 2 * beam_i_r2 * i_avg

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


# Generate bunch particles
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

# Create plasma stage
# -------------------
# Simulation box and grid parameters.
r_max = w_0 * 4  # Maximum radial extent of the box
r_max_plasma = w_0 * 3  # Maximum radial extent of the plasma.
l_box = 200e-6  # Length of simulation box
xi_max = xi_0 + 45e-6  # Right edge of the box in the speed-of-light frame
xi_min = xi_max - l_box  # Left edge of the box in the speed-of-light frame
res_beam_r = 5.  # default: 5, resolution we want for the beam radius
dr = np.min([sx0,sy0]) / res_beam_r  # make the radial resolution depend on the RMS beam size to avoid artifacts
dz = 1 / kp / 40
nr = int(r_max / dr)
nz = int(l_box / dz)

# Create plasma stages.
plasma_plateau = PlasmaStage(
    length=l_plateau,
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
    dz_fields= 1 * l_box,
    ppc=2,
    # parabolic_coefficient=rel_delta_n_over_w2,  # deprecated: use density_profile instead
    laser_envelope_substeps=4,
    laser_envelope_nxi=nz * 4,
    max_gamma=25,
    field_diags=["rho", "E", "a"],
)

if __name__ == "__main__":

    # Track the beam through the plasma stage
    beamline = Beamline([ plasma_plateau ])
    beamline.track(bunch, opmd_diag=True, show_progress_bar=True, diag_dir='diags')

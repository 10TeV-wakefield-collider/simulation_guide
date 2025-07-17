import sys
import os 
import numpy as np
import h5py
import matplotlib.pyplot as plt
import shutil
import importlib
import json
import subprocess, os 
import picmi_qpad as picmi 


cst = picmi.constants
codename = picmi.codename
#################### simulation reference density ####################
n0 = 1e16 * 1e6 # number density in units of 1/m^3
wp = np.sqrt(cst.q_e**2 * n0/(cst.ep0 * cst.m_e))
kp = wp/cst.c


#################### grid params (in units of c/wp) ####################
zmin,zmax = -6.6/kp, 3/kp
rmin, rmax = 0, 9.395847/kp
nr, nz = 256, 512
n_modes = 1


#################### time step params (in units of 1/wp) ####################
dt = 10/wp
tmax = 20000.1/wp
ndump_diag = 10 # dump every 10 timsteps


#################### MPI_config ####################
cpu_split = [1,16] # QPAD cpu split



#################### drive beam params ####################
beam1_charge = 3.193e-9 # charge [C]
beam1_sigmas = [7.285e-6, 7.285e-6, 2.549e-5] # sigmas [m]
beam1_centroid_position = [0,0,0] # centroid at (0,0,0) 


#################### trailing beam (#2) params ####################
beam2_charge = 0.956e-9 # charge [C]
beam2_sigmas = [7.285e-6, 7.285e-6, 1.274e-05] #sigmas [m]
beam2_centroid_position = [2e-6, 0, -3e-4] # offset 2 micron in x, behind driver by 300 um


#################### common beam params (shared by both beams) ####################
gamma = 19569.47 # Energy [10 GeV] 
beta= np.sqrt(2 * gamma)/kp # matched condition
rms_vel = [beam1_sigmas[0] * gamma/beta, beam1_sigmas[1] * gamma/beta, 0] # matched condition


#################### Plasma profile ####################
plasma_density = n0 # uniform density


#################### PICMI section ####################
solver_dict, beam_dist_dict, plasma_dist_dict, part_diag_dict= {}, {}, {}, {}
beam_layout_dict, plasma_layout_dict = {}, {}
field_diag_dict = {}
sim_dict = {codename + '_nodes' : cpu_split, codename + '_n0' : n0}

# layouts of beams+plasma
ppc_beam = [2, 1, 2] # (2,2) ppc along (r,z), second index skipped
ppc_plasma = [4, 1] # 4 part per cell in r, second index skipped
num_theta_beam = num_theta_plasma = 8   
	


#################### mkdir sim folder ####################
sim_dir= 'hosing_sims/' + codename + '_sim'
os.makedirs(sim_dir, exist_ok = True)


# initialize grid 
grid = picmi.CylindricalGrid(
			number_of_cells           = [nr, nz],
			lower_bound               = [0. , zmin],
			upper_bound               = [rmax, zmax],
			lower_boundary_conditions = ['open', 'open'],
			upper_boundary_conditions = ['open', 'open'],
			n_azimuthal_modes = n_modes,
			moving_window_velocity    = [0,cst.c])


solver = picmi.ElectromagneticSolver( grid = grid , **solver_dict)


## add drive beam
drive_dist = picmi.GaussianBunchDistribution(
	n_physical_particles = abs(int(beam1_charge/cst.q_e)),
	rms_bunch_size       = beam1_sigmas,
	rms_velocity         = [cst.c * x for x in rms_vel],
	centroid_position    = beam1_centroid_position,
	centroid_velocity    = [0, 0, gamma*cst.c], **beam_dist_dict )

drive_beam = picmi.Species( name = 'drive', particle_type = 'electron', initial_distribution = drive_dist)


## add trailing beam
trailing_dist = picmi.GaussianBunchDistribution(
	n_physical_particles = abs(int(beam2_charge/cst.q_e)),
	rms_bunch_size       = beam2_sigmas,
	rms_velocity         = [cst.c * x for x in rms_vel],
	centroid_position    = beam2_centroid_position,
	centroid_velocity    = [0, 0, gamma*cst.c], **beam_dist_dict )

trailing_beam = picmi.Species(name = 'trailing', particle_type = 'electron', initial_distribution = trailing_dist)


# add plasma density
plasma_dist = picmi.UniformDistribution(density = plasma_density)

plasma = picmi.Species(particle_type = 'electron', 
						name = 'plasma',
						initial_distribution = plasma_dist)

# beam/plasma layouts 
beam_layout_dict[picmi.codename + '_num_theta'] = num_theta_beam 
beam_layout = picmi.GriddedLayout(
		grid = grid,
		n_macroparticle_per_cell = ppc_beam, 
		**beam_layout_dict)


plasma_layout_dict[picmi.codename + '_num_theta'] = num_theta_plasma 
plasma_layout = picmi.GriddedLayout(
			grid = grid,
			n_macroparticle_per_cell = ppc_plasma, 
			**plasma_layout_dict)


field_diag = picmi.FieldDiagnostic(data_list = ['E','rho'],
	                                   grid = grid,
	                                   period = ndump_diag,
	                                   **field_diag_dict)

part_diag = picmi.ParticleDiagnostic(period = 1,
                                     species = [drive_beam, trailing_beam],
                                     **part_diag_dict)


sim = picmi.Simulation(solver = solver, verbose = 1,\
	 time_step_size = dt, max_time =tmax,  **sim_dict)


sim.add_species(species = drive_beam, layout = beam_layout)
sim.add_species(species = trailing_beam, layout = beam_layout)
sim.add_species(species = plasma, layout = plasma_layout)
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)


sim.write_input_file(sim_dir + '/qpinput.json')



# PATH to executable (qpad.e,...)
path_to_exe = ''   # <---- specify path to executable here
assert path_to_exe, Exception('Need to specify path to ' + picmi.codename + ' executable')

# optionally run the code use srun or mpirun 
if(True):
	env = dict(os.environ)
	procs = np.prod(cpu_split)
	f = open(sim_dir + '/output.txt', "w")
	f2 = open(sim_dir + '/stderr.txt', "w")
	subprocess.run(["mpirun", "-n", str(procs), path_to_exe + "qpad.e"],stdout =f,stderr =f2, cwd = sim_dir, env=env)
	# subprocess.run(["srun", "-n", str(procs), "-c", str(2), "--cpu_bind=cores", path_to_exe + "qpad.e"],stdout =f,stderr =f2, cwd = sim_dir, env=env)
	f.close()
	f2.close()





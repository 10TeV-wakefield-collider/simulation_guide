# Copyright 2022-2023 Thamine Dalichaouch, Frank Tsung
# QPAD FACET extension of PICMI standard

import picmistandard
import numpy as np
import re
import math
import json 
from json import encoder
from itertools import cycle
import periodictable
from decimal import Decimal
import importlib

importlib.reload(picmistandard)
encoder.FLOAT_REPR = lambda o: format(o, '.4f')

codename = 'QPAD'
picmistandard.register_codename(codename)


def to_scientific_notation(value, nprec = 7):
	return float(f"{value:.{nprec}e}")


class constants:
	c = 299792458.
	ep0 = 8.8541878128e-12
	mu0 = 4 * np.pi * 1e-7
	q_e = 1.602176634e-19
	m_e = 9.1093837015e-31
	m_p = 1.67262192369e-27

picmistandard.register_constants(constants)
# Species Class
class Neutral(picmistandard.PICMI_Species):
	"""
	QPAD-Specific Parameters
	
	### Beam-specific parameters ####

	QPAD_beam_evolution: boolean, optional
		Toggles beam evolution

	QPAD_quiet_start: boolean, optional
		If turned on, a set of image particles will be added to suppress the statistic noise. 
	
	QPAD_np: integer(3), optional
		Number of beam particles distributed along each direction. The product is the total no of particles.
	
	### Plasma-specific parameters ####

	QPAD_ppc: integer(2), optional
		Number of macroparticles per cell in a xi-slice.


	"""
	# initialization 
	def init(self, kw):
		
		# self.charge = self.charge_state * constants.q_e
		self.charge = -constants.q_e
		self.mass = constants.m_e
		m = re.match(r'(?P<iso>#[\d+])*(?P<sym>[A-Za-z]+)', self.particle_type)
		element = periodictable.elements.symbol(m['sym'])
		if(m['iso'] is not None):
			element = element[m['iso'][1:]]
		if(self.charge_state is not None):
			assert self.charge_state <= element.number, Exception('%s charge state not valid'%self.particle_type)
			try:
				element = element.ion[self.charge_state]
			except ValueError:
				# Note that not all valid charge states are defined in elements,
				# so this value error can be ignored.
				pass
		self.element = element.number
		self.ion_max = kw.pop(codename + '_ion_max', self.element)

		# set profile type
		
		if(isinstance(self.initial_distribution, UniformDistribution)):
			self.profile_type = 'neutral'
			self.push_type = 'robust'
		elif(isinstance(self.initial_distribution, AnalyticDistribution)):
			self.profile_type = 'neutral'
			self.push_type = 'robust'
		elif(isinstance(self.initial_distribution, PiecewiseDistribution)):
			self.profile_type = 'neutral'
			self.push_type = 'robust'
		else:
			print('Warning: Only Uniform, Analytic, and Piecewise distributions are currently supported.')


	def normalize_units(self):
		# normalized charge, mass, density
		self.q = self.charge/constants.q_e
		self.m = self.mass/constants.m_e



	def fill_dict(self, keyvals, if_lasers):
		if(if_lasers):
			keyvals['push_type'] = self.push_type + '_pgc'
		else:
			keyvals['push_type'] = self.push_type
		keyvals['q'] = self.q
		keyvals['m'] = self.m
		keyvals['element'] = self.element
		keyvals['ion_max'] = self.ion_max
		q_scale = np.abs(self.charge/constants.q_e)
		if(not isinstance(self.initial_distribution, FileDistribution)):
			if(self.density_scale is not None):
				keyvals['density'] = self.initial_distribution.norm_density *self.density_scale * q_scale
			else:
				keyvals['density'] = self.initial_distribution.norm_density * q_scale
			keyvals['density'] = to_scientific_notation(keyvals['density'])
		self.initial_distribution.fill_dict(keyvals)

	def activate_field_ionization(self,model,product_species):
		return


# Species Class
class Species(picmistandard.PICMI_Species):
	"""
	QPAD-Specific Parameters
	
	### Beam-specific parameters ####

	QPAD_beam_evolution: boolean, optional
		Toggles beam evolution

	QPAD_quiet_start: boolean, optional
		If turned on, a set of image particles will be added to suppress the statistic noise. 
	
	QPAD_np: integer(3), optional
		Number of beam particles distributed along each direction. The product is the total no of particles.
	
	### Plasma-specific parameters ####

	QPAD_ppc: integer(2), optional
		Number of macroparticles per cell in a xi-slice.


	"""
	# initialization 
	def init(self, kw):
		part_types = {'electron': [-constants.q_e, constants.m_e] ,\
		'positron': [constants.q_e, constants.m_e],\
		'proton': [constants.q_e, constants.m_p],\
		'anti-proton' : [-constants.q_e, constants.m_p]}

		if(self.particle_type in part_types):
			if(self.charge is None): 
				self.charge = part_types[self.particle_type][0]
			if(self.mass is None): 
				self.mass = part_types[self.particle_type][1]
		else:
			self.charge = self.charge_state * constants.q_e
			m = re.match(r'(?P<iso>#[\d+])*(?P<sym>[A-Za-z]+)', self.particle_type)
			element = periodictable.elements.symbol(m['sym'])
			if(m['iso'] is not None):
				element = element[m['iso'][1:]]
			if(self.charge_state is not None):
				assert self.charge_state <= element.number, Exception('%s charge state not valid'%self.particle_type)
				try:
					element = element.ion[self.charge_state]
				except ValueError:
					# Note that not all valid charge states are defined in elements,
					# so this value error can be ignored.
					pass
			self.element = element
			if self.mass is None:
				self.mass = element.mass*periodictable.constants.atomic_mass_constant

		# print(self.charge)

		# Handle optional args for beams 
		self.beam_evolution = kw.pop(codename + '_beam_evolution', True)
		self.quiet_start = kw.pop(codename + '_quiet_start', True)
		

		# set profile type
		if(isinstance(self.initial_distribution, GaussianBunchDistribution)):
			self.profile_type = 'beam'
			self.geometry = 'cartesian'
			self.push_type = 'reduced'
		elif(isinstance(self.initial_distribution, FileDistribution)):
			self.profile_type = 'beam'
			self.push_type = 'reduced'
			self.geometry = 'cartesian'
		elif(isinstance(self.initial_distribution, UniformDistribution)):
			self.profile_type = 'species'
			self.push_type = 'robust'
		elif(isinstance(self.initial_distribution, AnalyticDistribution)):
			self.profile_type = 'species'
			self.push_type = 'robust'
		elif(isinstance(self.initial_distribution, PiecewiseDistribution)):
			self.profile_type = 'species'
			self.push_type = 'robust'
		else:
			print('Warning: Only Uniform and Gaussian distributions are currently supported.')


	def normalize_units(self):
		# normalized charge, mass, density
		self.q = self.charge/constants.q_e
		self.m = self.mass/constants.m_e



	def fill_dict(self, keyvals, if_lasers):
		if(self.profile_type == 'beam'):
			keyvals['evolution'] = self.beam_evolution
			keyvals['quiet_start'] = self.quiet_start
			keyvals['geometry'] = self.geometry
		else:
			if(if_lasers):
				keyvals['push_type'] = self.push_type + '_pgc'
			else:
				keyvals['push_type'] = self.push_type
		keyvals['q'] = self.q
		keyvals['m'] = self.m

		q_scale = np.abs(self.charge/constants.q_e)
		if(not isinstance(self.initial_distribution, FileDistribution)):
			if(self.density_scale is not None):
				keyvals['density'] = self.initial_distribution.norm_density *self.density_scale * q_scale
			else:
				keyvals['density'] = self.initial_distribution.norm_density * q_scale
			keyvals['density'] = to_scientific_notation(keyvals['density'])
		self.initial_distribution.fill_dict(keyvals)

	def activate_field_ionization(self,model,product_species):
		return


picmistandard.PICMI_MultiSpecies.Species_class = Species
class MultiSpecies(picmistandard.PICMI_MultiSpecies):
	def init(self, kw):
		return
		# for species in self.species_instances_list:
		# 	print(species.name)



class GaussianBunchDistribution(picmistandard.PICMI_GaussianBunchDistribution):
	"""
	QPAD-Specific Parameters
	
	### Plasma-specific parameters ####

	self.profile: integer
		Specifies profile-type of uniform plasma.
		profile = 0 (uniform plasma or piecewise linear function in z), 12 (piecewise linear along r and z)
		Profiles are multiplicative f(r,z) = f(r)  * f(z)

	self.s, self.r: float array
		Specifies longitudinal coordinates of piecewise profile in z=s and r.

	self.fs, self.fz: float array
		Species normalized densities along coordinates specified by self.fs and self.fz.

	"""
	def init(self, kw):
		self.profile = ['gaussian', 'gaussian', 'gaussian']
		self.piecewise_s = kw.pop(codename + '_piecewise_s', None)
		self.piecewise_fs = kw.pop(codename + '_piecewise_fs', None)
		self.if_piecewise = False
		if(self.piecewise_fs is not None and self.piecewise_s is not None):
			self.if_piecewise = True
			# print('piecewise fs/s', self.piecewise_fs)
			self.profile = ['gaussian', 'gaussian', 'piecewise-linear']
		self.alpha = kw.pop(codename + '_alpha', None)

	def normalize_units(self,species, density_norm):
		# get charge, peak density in unnormalized units
		part_charge = species.charge
		total_charge = part_charge * self.n_physical_particles
		if(species.density_scale is not None):
			total_charge *= species.density_scale
		
		peak_density = total_charge/(part_charge * np.prod(self.rms_bunch_size) * (2 * np.pi)**1.5)


		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c


		# QPAD takes the spot size in normalized units (k_pe sigma), uth is divergence (sigma_{gamma * beta}), and ufl is fluid velocity (gamma* beta) )

		# normalized spot sizes
		for i in range(3):
			self.rms_bunch_size[i] *= k_pe 
			self.centroid_position[i] *= k_pe 
			self.rms_velocity[i] /= constants.c 
			self.centroid_velocity[i] /= constants.c

		if(self.if_piecewise):
			self.piecewise_s = [i * k_pe for i in self.piecewise_s]
		self.gamma = self.centroid_velocity[2]
		self.centroid_position[2] *= -1
		# normalized charge, mass, density
		self.q = species.charge/constants.q_e
		self.m = species.mass/constants.m_e
		self.norm_density = peak_density/density_norm

		self.tot_charge = total_charge/(-constants.q_e * density_norm * k_pe**-3)

	def fill_dict(self,keyvals):
		keyvals['profile'] = self.profile
		keyvals['gamma'] = to_scientific_notation(self.gamma)
		# keyvals['gauss_center'] = [to_scientific_notation(i) for i in self.centroid_position]
		centroid_ = [0, 0, self.centroid_position[2]]
		keyvals['gauss_center'] = [to_scientific_notation(i) for i in centroid_]
		keyvals['perp_offset_x'] = [0, to_scientific_notation(self.centroid_position[0]), 0]
		keyvals['perp_offset_y'] = [0, to_scientific_notation(self.centroid_position[1]), 0]
		keyvals['push_type'] = 'reduced'
		if(self.alpha is not None):
			keyvals['alpha'] = self.alpha
		# keyvals['total_charge'] = self.tot_charge
		# QPAD coordinate in xi = ct-z

		if(self.if_piecewise):
			keyvals['piecewise_fx3'] = self.piecewise_fs
			keyvals['piecewise_x3'] = self.piecewise_s
			rms_size_ = [to_scientific_notation(i) for i in self.rms_bunch_size]
			rms_size_[2] = "none"
			keyvals['gauss_sigma'] = rms_size_

			for j in range(2):
				keyvals['range' + str(j+1)] = [to_scientific_notation(-4 * self.rms_bunch_size[j] +centroid_[j]),\
				 to_scientific_notation(4 * self.rms_bunch_size[j] + centroid_[j])]
			keyvals['range' + str(3)] = [to_scientific_notation(np.min(self.piecewise_s)), to_scientific_notation(np.max(self.piecewise_s))]
		else:
			keyvals['gauss_sigma'] = [to_scientific_notation(i) for i in self.rms_bunch_size]
			for j in range(3):
				keyvals['range' + str(j+1)] = [to_scientific_notation(-4 * self.rms_bunch_size[j] +centroid_[j]),\
				 to_scientific_notation(4 * self.rms_bunch_size[j] + centroid_[j])]
		
		# if(self.tot_charge is not None):
		# 	keyvals['total_charge'] = self.tot_charge
		# keyvals['gauss_sigma'] = [to_scientific_notation(i) for i in self.rms_bunch_size]

		keyvals['uth'] = [to_scientific_notation(i) for i in self.rms_velocity]
		

class FileDistribution(picmistandard.base._ClassWithInit):
	"""
	QPAD-Specific Parameters
	
	### Plasma-specific parameters ####

	self.profile: integer
		Specifies profile-type of uniform plasma.
		profile = 0 (uniform plasma or piecewise linear function in z), 12 (piecewise linear along r and z)
		Profiles are multiplicative f(r,z) = f(r)  * f(z)

	self.z, self.r: array
		Specifies longitudinal coordinates of piecewise profile in z and r.

	self.fz, self.fr: array
		Species normalized densities along coordinates specified by self.fz and self.fr.

	QPAD_r_min, QPAD_r_max: float, optional
		Radial range (i.e. QPAD_r_min <= r <= QPAD_r_max) for particles in UniformDistribution. Only required when specifying transverse lower_bounds or upper_bounds.
 
	"""
	def __init__(self, filename = None, beam_center = [0, 0 ,0], file_center = [0, 0, 0], has_spin = False, **kw):
		self.filename = filename
		self.beam_center = beam_center
		self.file_center = file_center
		self.has_spin = has_spin
		self.npmax = kw.pop(codename + '_npmax', 2*10**6)
		self.handle_init(kw)
		


	def normalize_units(self,species, density_norm):
		# normalized charge, mass, density
		# self.q = species.charge/constants.q_e
		# self.m = species.mass/constants.m_e

		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c

		self.file_center= [k_pe * i for i in self.file_center]
		self.beam_center= [k_pe * i for i in self.beam_center]



	def fill_dict(self,keyvals):
		keyvals['filename'] = self.filename
		keyvals['beam_center'] = self.beam_center
		keyvals['file_center'] = self.file_center
		keyvals['push_type'] = 'boris'
		keyvals['npmax'] = self.npmax

class UniformDistribution(picmistandard.PICMI_UniformDistribution):
	"""
	QPAD-Specific Parameters
	
	### Plasma-specific parameters ####

	self.profile: integer
		Specifies profile-type of uniform plasma.
		profile = 0 (uniform plasma or piecewise linear function in z), 12 (piecewise linear along r and z)
		Profiles are multiplicative f(r,z) = f(r)  * f(z)

	self.z, self.r: array
		Specifies longitudinal coordinates of piecewise profile in z and r.

	self.fz, self.fr: array
		Species normalized densities along coordinates specified by self.fz and self.fr.

	QPAD_r_min, QPAD_r_max: float, optional
		Radial range (i.e. QPAD_r_min <= r <= QPAD_r_max) for particles in UniformDistribution. Only required when specifying transverse lower_bounds or upper_bounds.
 
	"""
	def init(self,kw):
		# default profile for uniform plasmas
		self.profile = ['uniform', 'uniform']
		


	def normalize_units(self,species, density_norm):
		# normalize plasma density
		self.norm_density = self.density/density_norm

		# normalized charge, mass, density
		# self.q = species.charge/constants.q_e
		# self.m = species.mass/constants.m_e

		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c

		# self.density_expression =  str(self.norm_density)
		# self.norm_density = 1.0

		for i in range(3):
			if(self.lower_bound[i] is not None):
				self.lower_bound[i] *= k_pe
			if(self.upper_bound[i] is not None):
				self.upper_bound[i] *= k_pe

		for i in range(3):
			self.rms_velocity[i] /= constants.c 

		# if(np.any(self.directed_velocity != 0.0)):
		# 	print('Warning: ' + codename + ' does not support directed velocity for Analytic Distributions.')

	def fill_dict(self,keyvals):
		back_str,front_str = construct_bounds(self.lower_bound,self.upper_bound)
		keyvals['profile'] = self.profile
		keyvals['uth'] = [to_scientific_notation(i) for i in self.rms_velocity]
		keyvals['density'] = to_scientific_notation(self.norm_density)


class PiecewiseDistribution(picmistandard.base._ClassWithInit):
	"""
	QPAD-Specific Parameters
	
	### Plasma-specific parameters ####

	self.profile: integer
		Specifies profile-type of uniform plasma.
		profile = 0 (uniform plasma or piecewise linear function in z), 12 (piecewise linear along r and z)
		Profiles are multiplicative f(r,z) = f(r)  * f(z)

	

	"""
	def __init__(self, density, lower_bound=[None, None, None], 
		upper_bound=[None, None, None], rms_velocity=[0.0, 0.0, 0.0], 
		directed_velocity=[0.0, 0.0, 0.0], fill_in=None, piecewise_s = [0.0], piecewise_fs = [1.0], **kw):
		
		self.density = density
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound
		self.rms_velocity = rms_velocity
		self.directed_velocity = directed_velocity
		self.fill_in = fill_in
		self.piecewise_s = piecewise_s
		self.piecewise_fs = piecewise_fs
		self.profile = ['uniform', 'piecewise-linear']
		self.handle_init(kw)
		


	def normalize_units(self,species, density_norm):
		# normalize plasma density
		self.norm_density = self.density/density_norm

		# normalized charge, mass, density
		# self.q = species.charge/constants.q_e
		# self.m = species.mass/constants.m_e

		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c

		self.piecewise_s = [k_pe * i for i in self.piecewise_s]
		self.piecewise_fs = [i/density_norm for i in self.piecewise_fs]
		for i in range(3):
			if(self.lower_bound[i] is not None):
				self.lower_bound[i] *= k_pe
			if(self.upper_bound[i] is not None):
				self.upper_bound[i] *= k_pe

		for i in range(3):
			self.rms_velocity[i] /= constants.c 

		# if(np.any(self.directed_velocity != 0.0)):
		# 	print('Warning: ' + codename + ' does not support directed velocity for Analytic Distributions.')


	def fill_dict(self,keyvals):
		keyvals['profile'] = self.profile
		keyvals['uth'] = [to_scientific_notation(i) for i in self.rms_velocity]
		keyvals['density'] = to_scientific_notation(self.norm_density)
		keyvals['piecewise_s'] = [to_scientific_notation(i) for i in self.piecewise_s]
		keyvals['piecewise_fs'] = [to_scientific_notation(i) for i in self.piecewise_fs]





class AnalyticDistribution(picmistandard.PICMI_AnalyticDistribution):
	"""
	QPAD-Specific Parameters
	
	### Plasma-specific parameters ####

	self.profile: integer
		Specifies profile-type of uniform plasma.
		profile = 13 (analytic functions x, y, z)
		Profiles are multiplicative f(r,z) = f(r)  * f(z)

	"""
	def init(self,kw):
		# default profile for uniform plasmas
		self.profile = ['analytic', 'analytic']
		if(np.any(self.momentum_expressions == None)):
			print('Warning: QPAD does not support momentum expressions for Analytic Distributions.')


	def normalize_units(self,species, density_norm):

		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c

		for i in range(3):
			if(self.lower_bound[i] is not None):
				self.lower_bound[i] *= k_pe
			if(self.upper_bound[i] is not None):
				self.upper_bound[i] *= k_pe

		self.density_expression = normalize_math_func(self.density_expression, density_norm)
		self.density_expression =  self.density_expression + '/' + str(density_norm)
		self.norm_density = 1.0

		for i in range(3):
			self.rms_velocity[i] /= constants.c 

		# if(np.any(self.directed_velocity != 0.0)):
		# 	print('Warning: ' + codename + ' does not support directed velocity for Analytic Distributions.')


	def fill_dict(self,keyvals):
		# if(self.lower_bound is not None)
		back_str,front_str = construct_bounds(self.lower_bound,self.upper_bound)
		keyvals['profile'] = self.profile
		keyvals['uth'] = [to_scientific_notation(i) for i in self.rms_velocity]
		# keyvals['math_func'] = front_str + self.density_expression + back_str
		keyvals['math_func'] =  self.density_expression


class ParticleListDistribution(picmistandard.PICMI_ParticleListDistribution):
	def init(self,kw):
		raise Exception('Particle list distributions not yet supported in QPAD')


# constant, analytic, or mirror fields not yet supported in QPAD
class ConstantAppliedField(picmistandard.PICMI_ConstantAppliedField):
	def init(self,kw):
		raise Exception("Constant applied fields are not yet supported in QPAD")

class AnalyticAppliedField(picmistandard.PICMI_AnalyticAppliedField):
	def init(self,kw):
		raise Exception("Analytic applied fields are not yet supported in QPAD")

class Mirror(picmistandard.PICMI_Mirror):
	def init(self,kw):
		raise Exception("Mirrors are not yet supported in QPAD")


class ElectromagneticSolver(picmistandard.PICMI_ElectromagneticSolver):
	"""
	QPAD-Specific Parameters
	
	QPAD_maximum_iterations: integer
		Number of iterations for predictor corrector solver.
	"""
	def init(self, kw):
		self.maximum_iterations = kw.pop(codename + '_maximum_iterations', None)
		if(self.maximum_iterations == None):
			print('Defaulting to n_iterations = 10 for predictor corrector')
			self.maximum_iterations = 10
	def fill_dict(self,keyvals):
		keyvals['iter_max'] = self.maximum_iterations
		keyvals['iter_reltol'] = 1e-3
		keyvals['iter_abstol'] = 1e-3
		keyvals['relax_fac'] = to_scientific_notation(1e-3 * (self.grid.dr/0.02)**2)

		
		
class ElectrostaticSolver(picmistandard.PICMI_ElectrostaticSolver):
	def init(self, kw):
		raise Exception('This feature is not supported. Please use the Electromagnetic solver.')

## Throw Errors if trying to use 1D/2D/3D cartesian grids with QPAD
class Cartesian1DGrid(picmistandard.PICMI_Cartesian1DGrid):
	def init(self, kw):
		raise Exception(codename + ' does not support this feature. Please specify a Cylindrical Grid.')

class Cartesian2DGrid(picmistandard.PICMI_Cartesian2DGrid):
	def init(self, kw):
		raise Exception(codename + ' does not support this feature. Please specify a Cylindrical Grid.')

class Cartesian3DGrid(picmistandard.PICMI_Cartesian3DGrid):
	def init(self, kw):
		raise Exception(codename + ' does not support this feature. Please specify a Cylindrical Grid.')


class CylindricalGrid(picmistandard.PICMI_CylindricalGrid):
	def init(self, kw):
		dims = 2

		# second check to make sure window moving forward at c (window speed doesn't actually matter for QPAD)

		# check for open boundaries at r_max
		if(self.upper_boundary_conditions[0] != 'open' or self.lower_boundary_conditions[1] != 'open' or self.upper_boundary_conditions[1] !='open'): 
			print('QPAD Defaulting to open boundaries in r and z-directions.')

		self.dr = np.abs(self.upper_bound[0]- self.lower_bound[0])/self.number_of_cells[0]
		self.dz = np.abs(self.upper_bound[1]- self.lower_bound[1])/self.number_of_cells[1]
		self.boundary = 'open'
		self.r = [self.lower_bound[0], self.upper_bound[0]]
		self.z = [self.lower_bound[1], self.upper_bound[1]]

	def power_of_two_check(self,n):
		return (n & (n-1) == 0) and n != 0

	def normalize_units(self, density_norm):
		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c

		#normalize coordinates 
		for i in range(2):
			self.r[i] *= k_pe
			self.z[i] *= k_pe
		self.dr *= k_pe
		self.dz *= k_pe

	def fill_dict(self,keyvals):
		keyvals['grid'] = self.number_of_cells
		keyvals['max_mode'] = self.n_azimuthal_modes
		box = {}
		# box['r'] = self.r
		box['r'] = [to_scientific_notation(i) for i in self.r]
		# QPAD 3D is in xi not z (multiply by -1 + reverse z coordinate)
		box['z'] = [to_scientific_notation(i) for i in [-self.z[1], -self.z[0]]]
		keyvals['box'] = box
		keyvals['field_boundary'] = self.boundary


class FileLayout(picmistandard.base._ClassWithInit):
	"""
	QPAD-Specific Parameters
	
	QPAD_np_per_dimension: integer array, optional
		Part per dim in each direction. (for beams only)

	QPAD_npmax: integer, optional
		Particle buffer size per MPI partition.

	QPAD_num_theta: integer, optional
		Number of particles in azimuthal direction. Defaults to 8 * n_azimuthal_modes.
	"""

	def __init__(self, grid = None, **kw):
		# n_macroparticles is required.
		self.profile_type = 'file'


	def fill_dict(self, keyvals,profile_type):
		keyvals['profile_type'] = self.profile_type


class PseudoRandomLayout(picmistandard.PICMI_PseudoRandomLayout):
	"""
	QPAD-Specific Parameters
	
	QPAD_np_per_dimension: integer array, optional
		Part per dim in each direction. (for beams only)

	QPAD_npmax: integer, optional
		Particle buffer size per MPI partition.

	QPAD_num_theta: integer, optional
		Number of particles in azimuthal direction. Defaults to 8 * n_azimuthal_modes.
	"""

	def init(self,kw):
		# n_macroparticles is required.
		assert self.n_macroparticles is not None, Exception('n_macroparticles must be specified when using PseudoRandomLayout with QPAD')
		self.profile_type = 'random'


	def fill_dict(self, keyvals,profile_type):
		keyvals['npmax'] = self.n_macroparticles * 2 
		keyvals['total_num'] = self.n_macroparticles
		keyvals['profile_type'] = self.profile_type
		keyvals['random_theta'] = False
		if(profile_type == 'beam'):
			if(self.n_macroparticles_per_cell is not None):
				keyvals['ppc'] = self.n_macroparticles_per_cell
		elif(profile_type == 'species'):
			raise Exception('PseudoRandomLayout not compatible with non-beam species')
				


class GriddedLayout(picmistandard.PICMI_GriddedLayout):
	"""
	QPAD-Specific Parameters

	QPAD_npmax: integer, optional
		Particle buffer size per MPI partition.

	QPAD_num_theta: integer, optional
		Number of particles in azimuthal direction. Defaults to 8 * n_azimuthal_modes.
	"""
	def init(self,kw):
		self.npmax = kw.pop(codename + '_npmax', 2*10**6)
		# assert len(self.n_macroparticle_per_cell) !=2, print('Warning: '+ codename + ' only supports 2-dimensions for n_macroparticle_per_cell')
		# setting profile type to standard
		self.profile_type = 'standard'
		self.num_theta = kw.pop(codename + '_num_theta', 1)
		# if(self.num_theta * self.n_macroparticle_per_cell[1] < 8 * self.grid.n_azimuthal_modes):
		# 	self.num_theta = int((8 * self.grid.n_azimuthal_modes)/self.n_macroparticle_per_cell[1])
		# 	print('Warning: total azimthal ppc increased to ' + str(self.num_theta * self.n_macroparticle_per_cell[1]))

	def fill_dict(self,keyvals,profile_type):
		keyvals['npmax'] = self.npmax
		# keyvals['profile_type'] = self.profile_type
		if(profile_type == 'beam'):
			keyvals['ppc'] = self.n_macroparticle_per_cell
		elif(profile_type == 'species' or profile_type == 'neutral'):
			keyvals['ppc'] = self.n_macroparticle_per_cell[:2]
		keyvals['num_theta'] = self.num_theta
		keyvals['random_theta'] = False



class Simulation(picmistandard.PICMI_Simulation):
	"""
	QPAD-Specific Parameters

	QPAD_n0: float, optional
		Plasma density [m^3] to normalize units.

	QPAD_nodes: int(2), optional
		MPI-node configuration

	QPAD_interpolation: str, optional
		Interpolation order (linear for QPAD).

	QPAD_read_restart: boolean, optional
		Toggle to read from restart files.
	
	QPAD_restart_timestep: integer, optional
		Specifies timestep if read_restart = True.

	QPAD_random_seed: integer, optional
		No of seeds for pseudo-random numbers. Defaults to 10.

	QPAD_interpolation: str, optional
		Interpolation order (e.g. linear).

	QPAD_algorithm: str, optional
		Type of algorithm (standard, pgc, etc). Defaults to standard.

	QPAD_timings: bool, optional
		Toggle to report timings. Turned off by default.

	"""
	def init(self,kw):
		# set verbose default
		if(self.verbose is None):
			self.verbose = 0
		assert self.time_step_size is not None, Exception('QPAD requires a time step size for the 3D loop.')

		if(self.particle_shape not in  ['linear']):
			print('Warning: Defaulting to linear particle shapes.')
			self.particle_shape = 'linear'

		self.cpu_split = kw.pop(codename + '_nodes', [1, 1])

		### QPAD differentiates between beams, neutrals, and plasmas (species)
		self.if_beam = []

		# no of neutrals
		self.if_neutral = []

		# no of species
		self.if_species = []

		# check if normalized density is specified
		self.n0 = kw.pop(codename + '_n0', None)

		# set number of seeds for pseudo-random numbers
		self.random_seed = kw.pop(codename + '_random_seed', 10)

		# set algorithm type (default is standard)
		self.algorithm = kw.pop(codename + '_algorithm', 'standard')

		# set timings (default is true)
		self.if_timing = kw.pop(codename + '_timings', False)

		# normalize simulation time
		if(self.n0 is not None):
			self.normalize_simulation()

		# check to read from restart files
		self.read_restart = kw.pop(codename + '_read_restart', False)
		self.restart_timestep = kw.pop(codename + '_restart_timestep', -1)
		if(self.read_restart):
			assert self.restart_timestep != -1, Exception('Please specify ' + codename + '_restart_timestep')

		# check if dumping restart files
		self.dump_restart = kw.pop(codename + '_dump_restart', False)
		self.ndump_restart = kw.pop(codename + '_ndump_restart', -1)
		if(self.dump_restart):
			assert self.ndump_restart != -1, Exception('Please specify' + codename + '_ndump_restart')


	def normalize_simulation(self):
		w_pe = np.sqrt(constants.q_e**2.0 * self.n0/(constants.ep0 * constants.m_e) ) 
		if(self.max_time is not None):
			self.max_time *= w_pe
		self.time_step_size *= w_pe
		self.solver.grid.normalize_units(self.n0)
	
	def add_species(self, species, layout, initialize_self_field = None):
		if(isinstance(species, MultiSpecies)):
			for spec in species.species_instances_list:
				picmistandard.PICMI_Simulation.add_species( self, spec, layout,
										  initialize_self_field )

				# handle checks for beams
				self.if_beam.append(spec.profile_type == 'beam')
				self.if_neutral.append(spec.profile_type == 'neutral')
				self.if_species.append(spec.profile_type == 'species')
				spec.normalize_units()
			if(self.n0 is not None):
					species.initial_distribution.normalize_units(spec, self.n0)

		else:
			picmistandard.PICMI_Simulation.add_species( self, species, layout,
										  initialize_self_field )
			if(self.n0 is not None):
				species.initial_distribution.normalize_units(species, self.n0)
				species.normalize_units()
			# handle checks for beams
			self.if_beam.append(species.profile_type == 'beam')
			self.if_neutral.append(species.profile_type == 'neutral')
			self.if_species.append(species.profile_type == 'species')
	def add_laser(self, laser, injection_method):
		picmistandard.PICMI_Simulation.add_laser(self, laser, injection_method)
		if(injection_method is not None):
			print('Antenna is not supported in QPAD. Initializating laser in box at t=0.')
			laser.focal_position[2] -= injection_method.position[2]
		if(self.n0 is not None):
			laser.normalize_units(self.n0)


			


	def fill_dict(self, keyvals):
		if(self.max_time is None):
			self.max_time = self.max_steps * self.time_step_size

		# fill grid and mpi params
		keyvals['nodes'] = self.cpu_split
		self.solver.grid.fill_dict(keyvals)

		# fill simulation time and dt
		keyvals['time'] = to_scientific_notation(self.max_time)
		keyvals['dt'] = to_scientific_notation(self.time_step_size)
		keyvals['interpolation'] = self.particle_shape


		if(self.n0 is not None):
			keyvals['n0'] = to_scientific_notation(self.n0 * 1.e-6) # in density in cm^{-3}
		keyvals['nbeams'] = int(np.sum(self.if_beam))
		keyvals['nspecies'] = int(np.sum(self.if_species))
		keyvals['nneutrals'] = int(np.sum(self.if_neutral))
		keyvals['nlasers'] = len(self.lasers)
		self.solver.fill_dict(keyvals)
		keyvals['dump_restart'] = self.dump_restart
		if(self.dump_restart):
			keyvals['ndump_restart'] = self.ndump_restart
		keyvals['read_restart'] = self.read_restart
		if(self.read_restart):
			keyvals['restart_timestep'] = self.restart_timestep
		keyvals['verbose'] = self.verbose
		keyvals['if_timing'] = self.if_timing
		keyvals['random_seed'] = self.random_seed
		keyvals['algorithm'] = self.algorithm

	def write_input_file(self,file_name):
		total_dict = {}

		# simulation object handled
		sim_dict = {}
		self.fill_dict(sim_dict)

		# beam objects 
		beam_dicts = []

		# species objects
		species_dicts = [] 

		# neutral objects
		neutral_dicts = []

		# lasers objects
		laser_dicts = []

		# field object
		field_dict = {}
		# iterate over species handle beams first
		for i in range(len(self.species)):
			spec = self.species[i]
			temp_dict = {}
			self.layouts[i].fill_dict(temp_dict,spec.profile_type)
			self.species[i].fill_dict(temp_dict, len(self.lasers) > 0)

			# fill in source term diagnostics
			diags_srcs = []
			for j in range(len(self.diagnostics)):
				diag = self.diagnostics[j]
				if(isinstance(diag,ParticleDiagnostic) and spec not in diag.species):
					continue
				temp_dict2 = {}
				self.diagnostics[j].fill_dict_src(temp_dict2)
				diags_srcs.append(temp_dict2)
			temp_dict['diag'] = diags_srcs
			if(self.if_beam[i]):
				beam_dicts.append(temp_dict)
			elif(self.if_neutral[i]):
				neutral_dicts.append(temp_dict)
			else:
				species_dicts.append(temp_dict)

		for i in range(len(self.lasers)):
			laser = self.lasers[i]
			temp_dict = {}
			self.lasers[i].fill_dict(temp_dict)
			laser_dicts.append(temp_dict)
			self.lasers[i].fill_dict_fld(temp_dict,self.diagnostics)

		diags_flds = []
		for i in range(len(self.diagnostics)):
			diag = self.diagnostics[i]
			temp_dict = {}
			if(isinstance(diag,ParticleDiagnostic)):
				continue
			self.diagnostics[i].fill_dict_fld(temp_dict)
			diags_flds.append(temp_dict)

		field_dict['diag'] = diags_flds

		total_dict['simulation'] = sim_dict
		if(len(beam_dicts) > 0):
			total_dict['beam'] = beam_dicts
		if(len(species_dicts) > 0):
			total_dict['species'] = species_dicts
		if(len(neutral_dicts) > 0):
			total_dict['neutrals'] = neutral_dicts
		if(len(laser_dicts) > 0):
			total_dict['laser'] = laser_dicts
		total_dict['field'] = field_dict
		with open(file_name, 'w') as file:
			json.dump(total_dict, file, indent =4)

	def step(self, nsteps = 1):
		raise Exception('The simulation step feature is not yet supported for ' + codename + '. Please call write_input_file() to construct the input deck.')

class FieldDiagnostic(picmistandard.PICMI_FieldDiagnostic):
	"""
	QPAD-Specific Parameters

	"""
	def init(self,kw):
		assert self.write_dir != '.', Exception("Write directory feature not yet supported.")
		assert self.period > 0, Exception("Diagnostic period is not valid")
		self.field_list = []
		self.source_list = []
		if('E' in self.data_list):
			self.field_list += ['er_cyl_m','ephi_cyl_m','ez_cyl_m']
		if('B' in self.data_list):
			self.field_list += ['br_cyl_m','bphi_cyl_m','bz_cyl_m']
		if('rho' in self.data_list):
			self.source_list += ['charge_cyl_m']
		if('J' in self.data_list):
			self.source_list += ['jr_cyl_m','jphi_cyl_m','jz_cyl_m']


		if('Ex' in self.data_list or 'Er' in self.data_list):
			self.field_list.append('er_cyl_m')
		if('Ey' in self.data_list or 'Ephi' in self.data_list):
			self.field_list.append('ephi_cyl_m')
		if('Ez' in self.data_list or 'Ez' in self.data_list):
			self.field_list.append('ez_cyl_m')

		if('Bx' in self.data_list or 'Br' in self.data_list):
			self.field_list.append('br_cyl_m')
		if('By' in self.data_list or 'Bphi' in self.data_list):
			self.field_list.append('bphi_cyl_m')
		if('Bz' in self.data_list or 'Bz' in self.data_list):
			self.field_list.append('bz_cyl_m')

		if('Jx' in self.data_list or 'Jr' in self.data_list):
			self.source_list.append('jr_cyl_m')
		if('Jy' in self.data_list or 'Jphi' in self.data_list):
			self.source_list.append('jphi_cyl_m')
		if('Jz' in self.data_list or 'Jz' in self.data_list):
			self.source_list.append('jz_cyl_m')

		# need to add to PICMI standard
		if('psi' in self.data_list):
			self.field_list += ['psi_cyl_m']


	def fill_dict_fld(self,keyvals):
		keyvals['name'] = self.field_list
		keyvals['ndump'] = self.period

	def fill_dict_src(self,keyvals):
		keyvals['name'] = self.source_list
		keyvals['ndump'] = self.period


# QPAD does not support electrostatic and boosted frame diagnostic 
class ElectrostaticFieldDiagnostic(picmistandard.PICMI_ElectrostaticFieldDiagnostic):
	def init(self,kw):
		raise Exception("Electrostatic field diagnostic not supported in QPAD")

class LabFrameParticleDiagnostic(picmistandard.PICMI_LabFrameParticleDiagnostic):
	def init(self,kw):
		raise Exception("Boosted frame diagnostics not support in QPAD")

class LabFrameFieldDiagnostic(picmistandard.PICMI_LabFrameFieldDiagnostic):
	def init(self,kw):
		raise Exception("Boosted frame diagnostics not support in QPAD")


## to be implemented	
class ParticleDiagnostic(picmistandard.PICMI_ParticleDiagnostic):
	"""
	QPAD-Specific Parameters

	QPAD_sample: integer, optional
		Dumps every nth particle.
	"""
	def init(self,kw):
		# print(kw)
		assert self.write_dir != '.', Exception("Write directory feature not yet supported.")
		assert self.period > 0, Exception("Diagnostic period is not valid")
		print('Warning: Particle diagnostic reporting momentum, position and charge data')
		self.sample = kw.pop(codename + '_sample', 1) 

	def fill_dict_fld(self,keyvals):
		pass

	def fill_dict_src(self,keyvals):
		keyvals['name'] = ["raw"]
		keyvals['ndump'] = self.period
		keyvals['psample'] = self.sample


class GaussianLaser(picmistandard.PICMI_GaussianLaser):
	def init(self, kw):
		self.profile = ['gaussian', 'polynomial']
		self.iteration = 3

	def normalize_units(self, density_norm):
		# normalize quantities to plasma density and skin depths
		w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
		k_pe = w_pe/constants.c


		self.k0 = self.k0/k_pe
		self.waist = k_pe * self.waist
		self.duration = w_pe * self.duration 
		for i in range(3):
			self.focal_position[i] *= k_pe 
			self.centroid_position[i] *= k_pe


	def fill_dict(self,keyvals):
		keyvals['profile'] = self.profile
		keyvals['a0'] = self.a0
		keyvals['k0'] = self.k0
		keyvals['w0'] = self.waist
		keyvals['iteration'] = self.iteration
		keyvals['focal_distance'] = self.focal_position[2]
		keyvals['t_rise'] = self.duration * 1.5275
		keyvals['t_fall'] = self.duration * 1.5275
		keyvals['t_flat'] = 0
		keyvals['lon_center'] = -self.centroid_position[2]

	def fill_dict_fld(self, keyvals,diagnostics):
		tt = []
		for diag in diagnostics:
			temp_dict = {}
			if(diag.data_list is not None):
				if('E' in diag.data_list or 'Ex' in diag.data_list or 'Ey' in diag.data_list or 'Ez' in diag.data_list):
					temp_dict['name'] = ['a_cyl_m']
					temp_dict['ndump'] = diag.period
					tt.append(temp_dict)
					break

		keyvals['diag'] = tt

		
class LaserAntenna(picmistandard.PICMI_LaserAntenna):
	def init(self, kw):
		return

class BinomialSmoother(picmistandard.PICMI_BinomialSmoother):
	def init(self, kw):
		print("Warning: QPAD has no BinomialSmoother. Skipping feature.")

def normalize_math_func(math_func, density_norm):
	w_pe = np.sqrt(constants.q_e**2 * density_norm/(constants.ep0 * constants.m_e) ) 
	k_pe = w_pe/constants.c

	## handle overlap with funcs and variable names
	funcs = ['exp','max', 'tan', 'sqrt','not','int','rect','step']
	funcs2 = ['e1p','ma1', '1an', 'sqr1', 'no1', 'in1', 'rec1', 's1ep']
	math_func2 = math_func[:]
	for i in range(len(funcs)):
		key1, key2 = funcs[i], funcs2[i]
		math_func2 = '{}'.format(math_func2).replace(key1,key2)


	math_func2 = '{}'.format(math_func2).replace('x', '(' + format_decimal(1.0/k_pe) +'* x)')
	math_func2 = '{}'.format(math_func2).replace('y', '(' + format_decimal(1.0/k_pe) +'* y)')
	math_func2 = '{}'.format(math_func2).replace('z', '(' + format_decimal(1.0/k_pe) +'* z)')
	math_func2 = '{}'.format(math_func2).replace('t', '(' + format_decimal(1.0/w_pe) +'* t)')
	for i in range(len(funcs)):
		key1, key2 = funcs2[i], funcs[i]
		math_func2 = '{}'.format(math_func2).replace(key1,key2)
	
	## handle overlap with funcs and variable names
	return math_func2

def format_decimal(decimal):
	str_out= '%.6e' % Decimal(str(decimal))
	return str_out

def construct_bounds(lower_bound,upper_bound):
	front_str =''
	back_str = ''
	coords = ['x','y','z']
	for i in range(3):
		if(upper_bound[i] is not None):
			if(len(front_str) > 0):
				front_str = front_str + ' && ' + coords[i] + '<= (' + format_decimal(upper_bound[i]) + ') ' 
			else:
				front_str = front_str + coords[i] + '<= (' + format_decimal(upper_bound[i]) + ') ' 

		if(lower_bound[i] is not None):
			if(len(front_str) > 0):
				front_str = front_str + ' && ' + coords[i] + '>= (' + format_decimal(lower_bound[i]) + ') ' 
			else:
				front_str = front_str + coords[i] + '>= (' + format_decimal(lower_bound[i]) + ') ' 
			

	front_str = 'if(' + front_str + ','
	back_str = ', 0 )'
	return back_str,front_str

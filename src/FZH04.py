import numpy as np
import tools21cm as t2c
import matplotlib.pyplot as plt
import scipy

import cosmolopy.perturbation as cp
import cosmolopy.distance as cd
import cosmolopy.density as cden
import cosmolopy.reionization as cr
import cosmolopy.constants as cc
import cosmolopy.parameters as cparam

import cosmo_const

def set_cosmology(cosm, baryonic_effects=False):
	if type(cosm) == str: cosmo_const.set_params(dataset=cosm)
	cosmo = {}
	cosmo['omega_M_0'] = t2c.const.Omega0 #0.3
	cosmo['omega_lambda_0'] = t2c.const.lam #0.7
	cosmo['omega_k_0'] = 1-cosmo['omega_M_0']-cosmo['omega_lambda_0']#0.0
	cosmo['h'] = t2c.const.h #0.7
	cosmo['n'] = cosmo_const.n_s
	cosmo['sigma_8'] = cosmo_const.sigma_8
	cosmo['omega_b_0'] = t2c.const.OmegaB #0.0462 
	cosmo['omega_n_0'] = cosmo_const.omg_n_0 #0.0 
	cosmo['N_nu'] = cosmo_const.N_nu #0 
	cosmo['Y_He'] = t2c.const.abu_he_mass #0.24 
	cosmo['baryonic_effects'] = baryonic_effects
	return cosmo


def density_threshold(cosmo=None, z=None, sig2=None, T_min=1e4, M_min=False, zeta=30, baryonic_effects=False):
	cosmo = set_cosmology(cosmo, baryonic_effects=baryonic_effects)

	mass_min = cp.virial_mass(T_min, z, **cosmo) # Return in solar masses
	r_min = cp.mass_to_radius(mass_min, **cosmo) # Return in Mpc
	sigma_min = cp.sigma_r(r_min, 0., **cosmo)[0]

	K = scipy.special.erfinv(1-1./zeta)
	if sig2 is None: sig2 = np.linspace(0,25,100)
	sd = cp.sig_del(T_min, z, passed_min_mass=M_min, **cosmo) # Returns sigma_min, delta_c
	delta_c = sd[1]
	delta_x = delta_c - np.sqrt(2*K**2*(sigma_min**2-sig2))
	return delta_x

def linear_fit_density_threshold(cosmo=None, z=10., T_min=1e4, zeta=30):
	sig2 = np.linspace(0,25,100)
	delta_x = density_threshold(cosmo=cosmo, z=z, sig2=sig2, T_min=T_min, zeta=zeta)
	B1, B0 = np.polyfit(sig2, delta_x, 1)
	return B0, B1

def m_dndm()



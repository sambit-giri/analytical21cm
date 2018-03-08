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
import FZH04

def set_cosmology(cosmo, baryonic_effects=False):
	if type(cosmo) == str: cosmo_const.set_params(dataset=cosmo)
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

def plot_FZH_fig1(cosmo=None, z=None, T_min=1e4, M_min=False, zeta=30, baryonic_effects=False):
	"""Plot figure 1 from FZH (2004ApJ...613....1F) (no quantitative tests).
	"""
	#cosmo = set_cosmology(cosmo)

	#mass_min = cp.virial_mass(T_min, z, **cosmo) # Return in solar masses
	#r_min = cp.mass_to_radius(mass_min, **cosmo) # Return in Mpc
	#sigma_min = cp.sigma_r(r_min, 0., **cosmo)[0]

	#K = scipy.special.erfinv(1-1./zeta)
	#sig2 = np.linspace(0,25,100)
	#sd = cp.sig_del(T_min, z, passed_min_mass=M_min, **cosmo) # Returns sigma_min, delta_c
	#delta_c = sd[1]
	#delta_x = delta_c - np.sqrt(2*K**2*(sigma_min**2-sig2))

	##B1, B0 = np.polyfit(sig2, delta_x, 1)
	#delta_x_fit = B0 + B1*sig2

	B0, B1  = FZH04.linear_fit_density_threshold(cosmo=cosmo, z=z, T_min=T_min, zeta=zeta)
	delta_x = FZH04.density_threshold(cosmo=cosmo, z=z, T_min=T_min, zeta=zeta)
	delta_x_fit = B0 + B1*sig2

	plt.plot(sig2, delta_x, c='k')
	plt.plot(sig2, delta_x_fit, '--', c='k')
	plt.xlabel('$\sigma^2$(m)')
	plt.ylabel('$\delta_x$')

	
	
def plot_FZH_fig2(cosmo=None, z=None, T_min=1e4, M_min=False, zeta=30, baryonic_effects=False):
	"""Plot figure 2 from FZH (2004ApJ...613....1F) (no quantitative tests).
	"""
	cosmo = set_cosmology(cosmo)
	rho_crit, rho_0 = cden.cosmo_densities(**cosmology)
	sd = cp.sig_del(T_min, z, passed_min_mass=M_min, **cosmo) #Returns sigma_min, delta_c
	cf = cp.collapse_fraction(*sd)

def plot_FZH_fig3(cosmo=None, z=None, T_min=1e4, zeta=np.array([[500.],[40.],[12.]]), baryonic_effects=False):
	"""Plot figure 3 from FZH (2004ApJ...613....1F) (no quantitative tests).
	"""
	cosmo = set_cosmology(cosmo)
	
	if z is None:
		dz = 0.1
		z = np.arange(20., 5. - 1.5 * dz, -1. * dz)

	# calculate ionized fraction from collapse fraction
	x_fcol = cr.ionization_from_collapse(z, zeta, T_min, **cosmo)
	#linestyle = ['-', ':', '--']
	color = ['r', 'g', 'b']
	plt.figure()
	plt.subplot(2,1,1)
	for i in range(len(color)):
		plt.semilogy(z, x_fcol[i], ls='--', color=color[i])
	plt.axhline(y=0.75)
	plt.xlim(5,20)
	plt.ylim(1e-4, 1)
	plt.title("Compare to figure 3 from FZH (2004ApJ...613....1F)")

	plt.subplot(2,1,2)
	for i in range(len(color)):
		plt.plot(z, x_fcol[i], ls='--', color=color[i])
	plt.axhline(y=0.75)
	plt.xlim(5,20)
	plt.ylim(0, 1)


if __name__ == "__main__":
	plot_FZH_fig1(z=12, zeta=30)
	plt.show()

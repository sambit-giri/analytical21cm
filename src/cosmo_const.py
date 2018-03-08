import numpy as np
import tools21cm as t2c

n_s     = 0.96
sigma_8 = 0.9
omg_n_0 = 0.0
N_nu    = 0.0

def set_new_params_in_tools21cm(h=None, OMm=None, OMl=None, OMb=None, Y_He=None):
	if Y_He is not None: t2c.set_abundance_helium(Y_He)
	if h is not None: t2c.set_hubble_h(h)
	if OMm is not None: t2c.set_omega_matter(OMm)
	if OMl is not None: t2c.set_omega_lambda(OMl)
	if OMb is not None: t2c.set_omega_baryon(OMb)

def set_params(dataset='Planck'):
	assert cosmo.lower() in ['planck', 'wmap+bao+h0', 'wmap']
	global n_s, sigma_8
	if dataset.lower() == 'planck':
		# Setting parameters from Planck 2015 results. XIII
		set_new_params_in_tools21cm(h=0.68, OMm=0.31, OMl=0.69, OMb=0.0486, Y_He=None)
		n_s, sigma_8 = 0.968, 0.816
	elif dataset.lower() == 'wmap':
		# Setting parameters from Komatsu+2011
		set_new_params_in_tools21cm(h=0.7, OMm=0.274, OMl=0.726, OMb=0.0455, Y_He=None)
		n_s, sigma_8 = 0.967, 0.811
	elif dataset.lower() == 'wmap+bao+h0':
		# Setting parameters from Komatsu+2011
		set_new_params_in_tools21cm(h=0.7, OMm=0.276, OMl=0.724, OMb=0.0458, Y_He=None)
		n_s, sigma_8 = 0.968, 0.816
	else: print "Unknown dataset name."
	print "New parameters set from", dataset, "results."




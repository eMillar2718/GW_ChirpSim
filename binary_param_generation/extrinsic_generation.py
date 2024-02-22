import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def redshift_to_luminosity(redshift):

    cosmo = FlatLambdaCDM(H0=67, Om0=0.3)
    luminosity_distance = cosmo.luminosity_distance(redshift)

    return luminosity_distance

def extrinsic_parameter_sampling(number_of_samples):

    extrinsic_population = {}
    extrinsic_populations = [] 

    for i in range(number_of_samples):

        uniform_sample = np.random.uniform(-1,1)


        phi_12 = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phi_12"] = phi_12

        phi_jl = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phi_jl"] = phi_jl 

        cosine_sample = np.arccos(uniform_sample) - np.pi/2
        extrinsic_population["dec"] = cosine_sample 

        sine_sample = np.arcsin(uniform_sample)
        sine_sample_shifted = sine_sample + np.pi / 2 
        extrinsic_population["theta_jn"] = sine_sample_shifted 

        psi = np.random.uniform(0, 2*np.pi)
        extrinsic_population["psi"] = psi 

        phase = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phase"] = phase 
        
        ra = np.random.uniform(0, 2*np.pi)
        extrinsic_population["ra"] = ra 

        extrinsic_populations.append(extrinsic_population)

    return extrinsic_populations

populations = extrinsic_parameter_sampling(20)


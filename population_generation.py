import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def intrinsic_parameter_sampling():

    pass
    #return mass_1, mass_2, a_1, a_2, tilt_1, tilt_2, redshift

def redshift_to_luminosity(redshift):

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    luminosity_distance = cosmo.luminosity_distance(redshift)

    return luminosity_distance

def extrinsic_parameter_sampling(number_of_samples):
     
    phi_12_values = []
    phi_jl_values = []
    theta_jn_values = []
    psi_values = []
    phase_values = []
    ra_values = []
    dec_values = []

    for i in range(number_of_samples):

        uniform_sample = np.random.uniform(-1,1)


        phi_12 = np.random.uniform(0, 2*np.pi)
        phi_12_values.append(phi_12)

        phi_jl = np.random.uniform(0, 2*np.pi)
        phi_jl_values.append(phi_jl)

        cosine_sample = np.arccos(uniform_sample) - np.pi/2
        dec_values.append(cosine_sample)

        sine_sample = np.arcsin(uniform_sample)
        sine_sample_shifted = sine_sample + np.pi / 2 
        theta_jn_values.append(sine_sample_shifted)

        psi = np.random.uniform(0, 2*np.pi)
        psi_values.append(psi)

        phase = np.random.uniform(0, 2*np.pi)
        phase_values.append(phase)
        
        ra = np.random.uniform(0, 2*np.pi)
        ra_values.append(ra)

     

    return phi_12_values, phi_jl_values, theta_jn_values, psi_values, phase_values, ra_values, dec_values 
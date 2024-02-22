import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import pandas as pd

def redshift_to_luminosity(redshift):

    cosmo = FlatLambdaCDM(H0=67, Om0=0.3)
    luminosity_distance = cosmo.luminosity_distance(redshift)

    return luminosity_distance.value

def extrinsic_parameter_sampling(number_of_samples):

    
    extrinsic_populations = [] 

    for i in range(number_of_samples):

        extrinsic_population = {}

        uniform_sample = np.random.uniform(-1,1)

        phi_12 = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phi_12"] = phi_12

        phi_jl = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phi_jl"] = phi_jl 


        sine_sample = np.arcsin(uniform_sample)
        sine_sample_shifted = sine_sample + np.pi / 2 
        extrinsic_population["theta_jn"] = sine_sample_shifted 

        psi = np.random.uniform(0, 2*np.pi)
        extrinsic_population["psi"] = psi 

        phase = np.random.uniform(0, 2*np.pi)
        extrinsic_population["phase"] = phase 
        
        ra = np.random.uniform(0, 2*np.pi)
        extrinsic_population["ra"] = ra 

        cosine_sample = np.arccos(uniform_sample) - np.pi/2
        extrinsic_population["dec"] = cosine_sample 

        extrinsic_populations.append(extrinsic_population)

    return extrinsic_populations

df_intrinsic = pd.read_hdf('/home/ethanmillar/University/GW_ChirpSim/binary_param_generation/all_intrins_samps.hdf5')

extrinsic_populations = extrinsic_parameter_sampling(len(df_intrinsic))

#print(df_intrinsic.loc[0:10])

df_extrinsic = pd.DataFrame(extrinsic_populations)

df_population = pd.concat([df_intrinsic, df_extrinsic], axis=1)

redshift = df_population["z"]

luminosity_distances = []

for i in range(len(redshift)):
    luminosity_distance = redshift_to_luminosity(redshift.loc[i])
    luminosity_distances.append(luminosity_distance)
                
df_population["luminosity_distance"] = luminosity_distances

df_population["m2"] = df_population["m1"] * df_population["q"] 

df_population.to_hdf('full_population_samples.hdf5', key='full_population_samples')


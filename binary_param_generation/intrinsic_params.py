import gwpopulation
import numpy as np

def get_p_m1(hyperpost_samp, n=1000):
    """
    Returns gwpopulation.models.mass.SinglePeakSmoothedMassDistribution from array of masses and mass ratios, with hyperposterior sample in 
    form of pandas data frame.
    
    Parameters
    ----------
    hyperpost_samp: dict
        1 hyperposterior sample, defines shape of mass population
        ['alpha', 'beta', 'mmax', 'mmin', 'lam', 'mpp', 'sigpp', 'delta_m']
    n: int
        number of points at which to evaluate distribution
    
    Returns
    -------
    p_m1: numpy array
        1D array of probability values given mass given the hyperposterior sample
        This is calculated from a Powerlaw plus Peak model with low mass smoothing.
    """

    alpha = hyperpost_samp['alpha']
    beta = hyperpost_samp['beta']
    mmin = hyperpost_samp['mmin']
    mmax = hyperpost_samp['mmax']
    lam = hyperpost_samp['lam']
    mpp = hyperpost_samp['mpp']
    sigpp = hyperpost_samp['sigpp']
    delta_m = hyperpost_samp['delta_m']

    masses= np.linspace(mmin,100.,n)
    qs= np.linspace(0.,1.,n)

    param_dict = {'mass_1':masses, 'mass_ratio':qs}


    mass_model = gwpopulation.models.mass.SinglePeakSmoothedMassDistribution()
    smoothing = mass_model.smoothing(masses, mmin=mmin, mmax=mmax, delta_m=delta_m)

    p_m1 = mass_model.p_m1(param_dict, **{'alpha':alpha, 'mmin':mmin, 'mmax':mmax, 'lam':lam, \
        'mpp':mpp, 'sigpp':sigpp})
    
    smooth_p_m1 = smoothing*p_m1
    normed_p_m1 =smooth_p_m1/np.trapz(smooth_p_m1,masses)

    return normed_p_m1, masses

def get_p_q(mass, hyperpost_samp, n=1000):
    """
    Parameters
    ----------
    mass: numpy array
        array of primary masses
    hyperpost_samp: dict
        1 hyperposterior sample, defines shape of mass population
        ['alpha', 'beta', 'mmax', 'mmin', 'lam', 'mpp', 'sigpp', 'delta_m']
    n: int
        number of points at which to evaluate distribution
    
    Returns
    -------
    p_q: numpy array
        1D array of probability values given q given the hyperposterior sample
        The mass ratio is given by a powerlaw mass ratio model with slope beta.
    """
    qs= np.linspace(0.,1.,n)
    param_dict = {'mass_1':mass, 'mass_ratio':qs}

    beta = hyperpost_samp['beta']
    mmin = hyperpost_samp['mmin']
    delta_m = hyperpost_samp['delta_m']

    q_model = gwpopulation.models.mass.SinglePeakSmoothedMassDistribution()
    p_q = q_model.p_q(param_dict, beta, mmin, delta_m)/np.trapz(q_model.p_q(param_dict, beta, mmin, delta_m), qs)

    return p_q/np.trapz(p_q, qs), qs

def get_p_z(hyperpost_samp, n=1000):
    """
    Parameters
    ----------
    z: numpy array
        array of redshifts
    hyperpost_samp: dict
        1 hyperposterior sample, defines shape of redshift population
        needs dict keys ['lamb']
    
    Returns
    -------
    p_z: numpy array
        1D array of probability values given z given the hyperposterior sample
        The mass ratio is given by a powerlaw mass ratio model with slope beta.
    """

    lamb = hyperpost_samp['lamb']
    z= np.linspace(0.,2.3,n)
    param_dict = {'redshift':z}

    z_model = gwpopulation.models.redshift.PowerLawRedshift()
    p_z = z_model.probability(param_dict, **{'lamb':lamb})

    return p_z/np.trapz(p_z, z), z

def get_p_chi(hyperpost_samp, n=1000):
    """
    Parameters
    ----------
    chi: numpy array
        array of component spin magnitudes
    hyperpost_samp: dict
        1 hyperposterior sample, defines shape of mass population
        ['mu_chi', 'sigma_chi', 'xi_spin', 'sigma_spin', 'lamb', 'amax']
    
    Returns
    -------
    p_chi: numpy array
        1D array of probability values given chi given the hyperposterior sample
        The mass ratio is given by a powerlaw mass ratio model with slope beta.
    """
    chi= np.linspace(0.,1.,n)

    amax = hyperpost_samp['amax']
    #conversion between hyperposterior sample values and gwpopulation/gwtc-3 pop paper values:
    alpha_chi, beta_chi, amax = gwpopulation.conversions.mu_var_max_to_alpha_beta_max(hyperpost_samp['mu_chi'],hyperpost_samp['sigma_chi'], amax)

    param_dict = {'a_1':chi, 'a_2':chi}
    p_chi = gwpopulation.models.spin.iid_spin_magnitude_beta(param_dict, amax=amax, alpha_chi=alpha_chi, beta_chi=beta_chi)

    return p_chi/np.trapz(p_chi, chi), chi

def get_p_costilt(hyperpost_samp, n=1000):
    """
    Parameters
    ----------
    cos_tilt: numpy array
        array of componant cosine spin tilts
    hyperpost_samp: dict
        1 hyperposterior sample, defines shape of mass population
        ['mu_chi', 'sigma_chi', 'xi_spin', 'sigma_spin', 'lamb', 'amax']
    
    Returns
    -------
    p_chi: numpy array
        1D array of probability values given chi given the hyperposterior sample
        The mass ratio is given by a powerlaw mass ratio model with slope beta.
    """
    cos_tilt = np.linspace(-1.,1.,n)

    xi_spin = hyperpost_samp['xi_spin']
    sigma_spin = hyperpost_samp['sigma_spin']

    param_dict = {'cos_tilt_1':cos_tilt, 'cos_tilt_2':cos_tilt}
    p_costilt = gwpopulation.models.spin.iid_spin_orientation_gaussian_isotropic(param_dict, xi_spin, sigma_spin)


    return p_costilt/np.trapz(p_costilt, cos_tilt), cos_tilt

def CDF(distr, theta):
    """
    Calculates a CDF of a PDF (distr) supplied at points theta
    """
    CDF = []
    for i in range(len(theta)):
        CDF.append(np.trapz(distr[:i+1], theta[:i+1]))
    return CDF, theta

def sample_1D(distr, theta, N):
    """
    returns samples from distribution calculated by interpolation of a given distribution

    Parameters
    ----------
    distr: numpy array
        values of a probability distribution of one of the binary parameters
    theta: numpy
        values of binary parameter at locations of probability distribution
    N: int
        number of samples to return from distr

    Returns
    -------
    samps: numpy array
        samples from distr
    """
    rand = np.random.random(N)
    CDF_theta, theta = CDF(distr, theta)
    samps = np.interp(rand, CDF_theta, theta)
    return samps

def sample_intrinsic(hyperpost_samp, N):

    sample = {'m1':0, 'q':0, 'z':0, 'chi_1':0, 'chi_2':0, 'costilt_1':0, 'costilt_2':0}

    #mass model
    p_m1, masses = get_p_m1(hyperpost_samp)
    sample['m1'] = sample_1D(p_m1, masses, N)


    #mass ratio
    for m1 in sample['m1']:
        p_q, qs = get_p_q(np.array([m1]), hyperpost_samp)
        q = sample_1D(p_q, qs, N)
    sample['q'] = q

    #redshift
    p_z, zs = get_p_z(hyperpost_samp)
    sample['z'] = sample_1D(p_z, zs, N)

    #spin magnitudes
    #should chi_1 and chi_2 be sampled from different spin populations or the same population?
    p_chi, chis = get_p_chi(hyperpost_samp)
    sample['chi_1'] = sample_1D(p_chi, chis, N)

    #spin tilts
    #should chi_1 and chi_2 be sampled from different spin populations or the same population?
    p_costilt, costilts = get_p_costilt(hyperpost_samp)
    sample['costilt_1'] = sample_1D(p_costilt, costilts, N)

    return sample

def sample_m1(hyperpost_samp, N):
    p_m1, masses = get_p_m1(hyperpost_samp)
    samples = sample_1D(p_m1, masses, N)
    return samples
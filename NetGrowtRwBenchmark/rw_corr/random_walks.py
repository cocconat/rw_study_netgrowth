import numpy as np
import sys

def correlated_rw(params,probability_function,_sigma=None):
    """
    Run a correlated random walk with a certain probability distribution

    Parameters:
    -----------------

    params['pers_length'] = 100
    params['seed']= np.random.randint(1000)
    params['exp_len']= 10000
    params['f'] = np.exp(-1. / params['pers_length'])
    params['f_sqrt'] = np.sqrt(1-f*f)
    params['sigma']= 0.143044444
    params['theta_0']= 0
    params['resolution']=101

    probability_function: (function) discrete or continuos probability.

    Returns:
    -----------
    np.array path: sequence of float.

    """
    f, f_sqrt, sigma, theta_0 = params['f'], params['f_sqrt'], params['sigma'],params['theta_0'],
    theta, r = theta_0,0
    if _sigma is not None:
        sigma =_sigma
        print (sigma)
    path=[]
    s=0
    arrays = generate_discrete_gaussian(params['resolution'])
    for _ in range(params['exp_len']):
        path.append(theta)
    ##update
        r = f*r + f_sqrt*s
        s =probability_function(arrays)
        theta = theta + sigma*(r+s)
    return np.array(path)

def correlated_gaussian(params,probability_function,_sigma=None):

    """
    Run a process of correlated gaussian with a certain probability distribution

    Parameters:
    -----------------

    params['pers_length'] = 100
    params['seed']= np.random.randint(1000)
    params['exp_len']= 10000
    params['f'] = np.exp(-1. / params['pers_length'])
    params['f_sqrt'] = np.sqrt(1-f*f)
    params['sigma']= 0.143044444
    params['theta_0']= 0
    params['resolution']=101

    probability_function: (function) discrete or continuos probability.

    Returns:
    -----------
    np.array path: sequence of float.

    """
    f, f_sqrt, sigma, theta_0 = params['f'], params['f_sqrt'], params['sigma'],params['theta_0'],
    theta, r = theta_0,0
    if _sigma is not None:
        sigma =_sigma
        print (sigma)
    path=[]
    arrays = generate_discrete_gaussian(params['resolution'])
    for _ in range(params['exp_len']):
        path.append(theta)
    ##update
        r = f*r + f_sqrt*probability_function(arrays)
        theta = theta_0 + sigma*r
    return np.array(path)

def classic_rw(params, probability_function,  _sigma=None):
    _, _, sigma,theta_0 = params['f'], params['f_sqrt'], params['sigma'],params['theta_0'],
    theta, r = theta_0,0
    if _sigma is not None:
        sigma =_sigma
        print (sigma)
    path=[]
    arrays = generate_discrete_gaussian(params['resolution'])
    for _ in range(params['exp_len']):
        path.append(theta)
    ##update
        r =probability_function(arrays)
        theta = theta + sigma*r
    return np.array(path)

def generate_discrete_gaussian(resolution):
    angles=[]
    space = 4./(resolution/2. -0.5)
    for x in range(resolution):
        angles.append(space*(x+(1 -resolution)/2))
    weights=[]
    for x in angles:
        weights.append( space/np.sqrt(2*3.14) * (0.5*np.exp(-0.5*(x-space/2.)**2.)+0.5*np.exp(-0.5*(x+space/2.)**2.)))
    print (sum(weights))
    return angles, weights/(sum(weights))


def initialize_exp(hardcoded =False):
    params={}
    if hardcoded == False:
        params['pers_length'] = int(sys.argv[1])
        params['seed']=int(sys.argv[2])
        params['exp_len']=int(sys.argv[3])
    else:
        params['pers_length'] = 100
        params['seed']= np.random.randint(1000)
        params['exp_len']= 10000
    params['f'] = np.exp(-1. / params['pers_length'])
    f = np.exp(-1. / params['pers_length'])
    params['f_sqrt'] = np.sqrt(1-f*f)
    params['sigma']= 0.143044444
    params['theta_0']= 0
    params['resolution']=101
    return params

def variance_changing(params):
    status = params['theta_0'], 0
    all_paths=[]
    for sigma in np.arange(0.2, 2.2, 0.3):
        path=[]
        np.random.seed(params['seed'])
        for _ in range(params['exp_len']):
            status =classic_rw(params,status, sigma)
            path.append(status)
        all_paths.append((np.array(path), str(sigma)))
    return all_paths

def run_walkers(params):
    """
    Run some algorithm with the defined parameters
    """
    print ( "correlation_coeff {}, sigma {}, theta_0 {}".format(params['f'],
                                                    params['sigma'],
                                                    params['theta_0']))
    np.random.seed(params['seed'])
    path_correlated=correlated_rw(params,continuos_probability)
    np.random.seed(params['seed'])
#     path_discrete=correlated_rw(params,discrete_probability)
    # np.random.seed(params['seed'])
    path_classic_rw=classic_rw(params,continuos_probability)

    return [(path_correlated,'correlated'), (path_classic_rw, 'uncorrelated')]
    #return [(path_classic_rw, 'random_walk')]

def discrete_probability(arrays):
    """
    return
    """
    delta = np.random.choice(arrays[0],p=arrays[1])
    return delta

def continuos_probability(arrays):
    return np.random.normal()



def normalize(array):
    norm = np.sum(array)
    return array/norm

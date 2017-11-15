from MeasurePersistence import Test
import numpy as np


def LinearSearch(neuron_params):
    msd=0
    cosine=0
    tortuosity=0
    max_attempt =100
    attempt = 0
    neuron_params['rw_memory_tau'] = 0.1
    old_cosine = 1
    while abs(cosine - min_persistence) > 20 or attempt > max_attempt:
        if cosine < min_persistence:
            neuron_params['rw_memory_tau'] = neuron_params['rw_memory_tau']*1.2
        elif cosine > min_persistence:
            neuron_params['rw_memory_tau'] = neuron_params['rw_memory_tau']*0.8
        attempt +=1
        if neuron_params['rw_memory_tau']< 0.1 or neuron_params['rw_memory_tau']>1000 or abs(old_cosine-cosine)<0.001:
            break
        old_cosine = cosine
        msd, tortuosity, cosine= Test(neuron_params)

    min_cosine = cosine
    min_msd = msd
    min_tortuosity = tortuosity
    min_memory = neuron_params['rw_memory_tau']

    cosine=0
    attempt = 0
    old_cosine = 1
    while abs(cosine - max_persistence) >20 or attempt > max_attempt:
        if cosine < max_persistence:
            neuron_params['rw_memory_tau'] = neuron_params['rw_memory_tau']*1.2
        elif cosine > max_persistence:
            neuron_params['rw_memory_tau'] = neuron_params['rw_memory_tau']*0.8
        attempt +=1
        if neuron_params['rw_memory_tau']< 0.1 or neuron_params['rw_memory_tau']>1000 or abs(old_cosine-cosine)<0.001:
            break
        old_cosine = cosine
        msd, tortuosity, cosine= Test(neuron_params)

    max_cosine = cosine
    max_msd = msd
    max_tortuosity = tortuosity
    max_memory = neuron_params['rw_memory_tau']

    return max_memory, min_memory, max_cosine, min_cosine, max_msd,min_msd, max_tortuosity, min_tortuosity

## Grid search on the memory
couple=[]
sigma_precision=20
corr_precision=20
neuron_params = {
        "axon_angle":0.,
        "use_tubulin": False,
        "rw_delta_corr": 0.001,
        "rw_memory_tau": 10.,
        "rw_sensing_angle":0.05,
        "speed_growth_cone": 1.,
        }

max_persistence = 300
min_persistence = 50

for g in np.linspace(0.1,2,11):
    for x in np.linspace(0.1,1,11):
        neuron_params['rw_sensing_angle']=x
        neuron_params['rw_delta_corr']=g
        try:
            values = LinearSearch(neuron_params)
            couple.append((x,g,*values))
        except:
            couple.append((x,g,*tuple([np.nan]*len(values))))
            pass


data_array = np.array(couple)
np.savetxt("plausible_values",data_array)
# data_array = data_array.reshape()

def pcplot(data_array):
    grid_corr = data_array[:,1].reshape(sigma_precision,corr_precision)

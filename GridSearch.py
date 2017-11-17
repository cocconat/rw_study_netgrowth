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
if __name__=="__main__":
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

def pcplot(data):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    fig, ((ax1, ax2),( ax3, ax4)) = plt.subplots(2,2)
    ax1.set_xlabel("sigma")
    ax1.set_ylabel("correlation coefficient")
    ax2.set_xlabel("sigma")
    ax2.set_ylabel("correlation coefficient")
    ax3.set_xlabel("sigma")
    ax3.set_ylabel("correlation coefficient")
    ax4.set_xlabel("sigma")
    ax4.set_ylabel("correlation coefficient")
    a=ax3.pcolor(data[:,:7,0],np.exp(1./-data[:,:7,1]),np.exp(1./-data[:,:7,2]))
    b=ax4.pcolor(data[:,:7,0],np.exp(1./-data[:,:7,1]),np.exp(1./-data[:,:7,3]))
    c=ax1.pcolor(data[:,:7,0],np.exp(-data[:,:7,1]),-data[:,:7,4])
    d=ax2.pcolor(data[:,:7,0],np.exp(-data[:,:7,1]),-data[:,:7,5])
    fig.colorbar(a)
    fig.colorbar(b)
    fig.colorbar(c)
    fig.colorbar(d)
    fig.tight_layout()

    ax1.set_title(r"$\langle b_i, b_j \rangle$ for max memory coeff.")
    ax2.set_title(r"$\langle b_i, b_j \rangle$ for min memory coeff.")
    ax3.set_title("Memory coefficient for 300 $\mu m$")
    ax4.set_title("Memory coefficient for 50 $\mu m$")


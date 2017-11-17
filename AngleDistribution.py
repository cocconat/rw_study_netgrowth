import NetGrowth
from NetGrowthRwBenchmark import AnalyseNetgrowthRW, CleanFolder
import numpy as np
import matplotlib.pyplot as plt
import os

def RunNetGrowth(n_samples, sim_length, n_procs, neuron_params, save_path = "tmp_measure", plot=False):
    """
    Run NetGrowth simulation
    """
    kernel={"seeds":[33,57,19,37,79,87,11][:n_procs],
            "num_local_threads":n_procs,
            "resolution":1.}
    experiment_params={}
    experiment_params["num_neurons"]=n_samples
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())
    culture_file =  "../culture/angle20.svg"
    culture = NetGrowth.CreateEnvironment(culture_file, min_x=0, max_x=1000)
    pos_left = culture.seed_neurons(neurons=experiment_params["num_neurons"], xmax=200, soma_radius=10.)

    neuron_params['growth_cone_model']='random_walk'
    neuron_params['position'] = pos_left

    gids =None
    # plt.show()
    gids = NetGrowth.CreateNeurons( experiment_params["num_neurons"],
                                        "random_walk",
                                        culture=culture,
                                        params=neuron_params,
                                        num_neurites=1
                                        )
    NetGrowth.Simulate(sim_length)
    fig, ax = plt.subplots()
    NetGrowth.plot.PlotNeuron(gid=range(experiment_params["num_neurons"]), culture=culture, soma_color="k",
                       axon_color='g', axis=ax, show=True)
    # if plot:
        # NetGrowth.PlotNeuron()
    NetGrowth.SaveJson(filepath=save_path)
    NetGrowth.SaveSwc (filepath=save_path,swc_resolution = 10)
    # NetGrowth.SaveJson(filepath=tmp_dir)
    NetGrowth.SaveSwc(filepath=os.path.join(os.getcwd(),save_path),swc_resolution = 10)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()


def Test(neuron_params, sim_length=500, sim_samples=30, plot=False):
    folder = os.path.join(os.getcwd(),"tmp_measure")
    CleanFolder(folder)
    RunNetGrowth(sim_samples, sim_length,5,  neuron_params, folder )
    NG_populations = NetGrowth.SimulationsFromFolder(folder)
    ensembles, _ =AnalyseNetgrowthRW(NG_populations)
    return ensembles
    # # return ensembles
    # # rw_corr.plot_results(ensembles, plot=True)
    # info =InfoFromJson(os.path.join(folder,"info.json"))
    # fit = fits[list(fits.keys())[0]]
    # fit = analyse_fit(fit, info=info)
    # fit= OnlyValues(fit)
    # if plot:
        # CleanFolder(folder)
        # RunNetGrowth(1, 1, neuron_params, folder,True)
    # CleanFolder(folder,make=False)
    # print(" ################################### \n")
    # print(" Memory Tau: {} um \n".format(fit["memory"]))
    # print(" Correlation Tau: {} um \n".format(fit["pers_gauss"]))
    # print(" Sigma: {} \n".format(fit["sigma"]))
    # print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    # print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    # print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    # print(" ################################## \n")
    # return fit["pers_length"], fit["tortuosity_local"], fit['cosine']

    # # json.dump(fit,open(os.path.join(folder,"fit.json"),'w'))


if __name__ =="__main__":
    neuron_params = {
        "filopodia_wall_affinity": 5.,
        "filopodia_finger_length": 20.,
        "filopodia_angular_resolution": 30,
        "axon_angle":0.,
        "use_tubulin": False,
        "rw_delta_corr": 0.1,
        "rw_memory_tau": 0.7,
        "rw_sensing_angle":0.15,
        "speed_growth_cone": 1.05,
        }
    ensemble1=Test(neuron_params)
    neuron_params["filopodia_wall_affinity"]=10.
    ensemble2=Test(neuron_params)
    neuron_params["filopodia_wall_affinity"]=1.
    ensemble3=Test(neuron_params)
    plt.hist(ensemble1[0].theta[:,-1])
    plt.hist(ensemble2[0].theta[:,-1], alpha=0.3)
    plt.hist(ensemble3[0].theta[:,-1], alpha=0.3)
    plt.show()


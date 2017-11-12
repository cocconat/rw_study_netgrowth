#!/usr/bin/env python
#-*- coding:utf-8 -*-

# MeasurePersistence measuresthe persistence length and the tortuosity resulting from a set of given parameters, this is very useful for NetGrowth user.


#These parameters are ok with a persistence length shoerter than 1mm
sim_length = 2000
max_len    = 200

#These parameters must be set.
neuron_params = {
        "axon_angle":0.,
        "use_tubulin": False,
        "rw_delta_corr": 20.,
        "rw_memory_tau": 500.,
        "rw_sensing_angle":0.029433,
        "speed_growth_cone": 9.95,
        }

# This software is part of NetGrowth project and SENEC initiative.


import NetGrowth
from  NetGrowthRwBenchmark import AnalyseNetgrowthRW, analyse_fit
import os, json, shutil
import numpy as np


def InfoFromJson(info_file, dump_kernel_params=False, dump_neuron_params=True):
    """
    Return neurons and kernel parameters into dictionary, from json file.

    Parameters:
    ==========

    info_file: info.json file produced by NetGrowth
    kernel_params: bool, retrieve Kernel parameters into

    """
    experiment_params = json.load(open(info_file, 'r',encoding='UTF8'))
    neuron_params = experiment_params['neurons']['0']['axon_params']
    kernel_params = experiment_params['kernel']
    if dump_kernel_params:
        return neuron_params, kernel_params
    else:
        return neuron_params



def RunNetGrowth(n_samples, n_procs, neuron_params, save_path = "tmp_measure", plot=False):
    kernel={"seeds":[33,57,19,37,79,87,11][:n_procs],
            "num_local_threads":n_procs,
            "environment_required":False}
    experiment_params={}
    experiment_params["num_neurons"]=n_samples
    initial_position=[0,0]
    positions=[]
    for x in range(experiment_params['num_neurons']):
        positions.append(initial_position)
    neuron_params["position"]   =positions
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())

    neuron_params['growth_cone_model']='random_walk'

    gids =None
    # plt.show()
    gids = NetGrowth.CreateNeurons(     experiment_params["num_neurons"],
                                        "random_walk",
                                        params=neuron_params,
                                        axon_params=neuron_params,
                                        dendrite_params=None,
                                        num_neurites=1,
                                        set_position=True)
    NetGrowth.Simulate(sim_length)
    if plot:
        NetGrowth.PlotNeuron()
    NetGrowth.SaveJson(filepath=save_path)
    NetGrowth.SaveSwc (filepath=save_path,swc_resolution = 10)
    # NetGrowth.SaveJson(filepath=tmp_dir)
    NetGrowth.SaveSwc(filepath=os.path.join(os.getcwd(),save_path),swc_resolution = 10)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()

def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)

def OnlyValues(main_dict):
    values={}
    for key in main_dict.keys():
        try:
            values[key]=main_dict[key]["values"]['a0']
        except:
            values[key]=main_dict[key]
    return values

def Test(neuron_params):
    folder = os.path.join(os.getcwd(),"tmp_measure")
    CleanFolder(folder)
    RunNetGrowth(100, 5, neuron_params, folder )
    NG_populations = NetGrowth.SimulationsFromFolder(folder)
    ensembles, fits =AnalyseNetgrowthRW(NG_populations,int(max_len))
    # return ensembles
    # rw_corr.plot_results(ensembles, plot=True)
    info =InfoFromJson(os.path.join(folder,"info.json"))
    fit = fits[list(fits.keys())[0]]
    fit = analyse_fit(fit, info=info)
    fit= OnlyValues(fit)
    # CleanFolder(folder)
    # RunNetGrowth(1, 1, neuron_params, folder,1)
    # CleanFolder(folder,make=False)
    print(" ################################### \n")
    print(" Memory Tau: {} um \n".format(fit["memory"]))
    print(" Correlation Tau: {} um \n".format(fit["pers_gauss"]))
    print(" Sigma: {} \n".format(fit["sigma"]))
    print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    print(" ################################## \n")
    return fit["pers_length"], fit["tortuosity_local"], fit

    # json.dump(fit,open(os.path.join(folder,"fit.json"),'w'))


Test(neuron_params)

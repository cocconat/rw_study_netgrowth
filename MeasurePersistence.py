#!/usr/bin/env python
# -*- coding:utf-8 -*-
# This software is part of the NetGrowth project and the SENEC initiative

# MeasurePersistence measures the persistence length and the tortuosity
# resulting from a set of given parameters, this is very useful for NetGrowth
# user.

import NetGrowth
from NetGrowthRwBenchmark import AnalyseNetgrowthRW, analyse_fit, OnlyValues
import os
import json
import shutil
import numpy as np

# These parameters are ok with a persistence length shoerter than 1mm
sim_length = 1000
max_len = 200

# These parameters must be set.
neuron_params = {
    "axon_angle": 0.,
    "use_tubulin": False,
    "rw_delta_corr": 0.1,
    "rw_memory_tau": 0.7,
    "rw_sensing_angle": 0.15,
    "speed_growth_cone": 1.05,
}


def InfoFromJson(info_file, dump_kernel_params=False, dump_neuron_params=True):

    """
    Return neurons and kernel parameters into dictionary, from json file.

    Parameters:
    ==========
    info_file: info.json file produced by NetGrowth
    kernel_params: bool, retrieve Kernel parameters into
    """
    experiment_params = json.load(open(info_file, 'r', encoding='UTF8'))
    neuron_params = experiment_params['neurons']['0']['axon_params']
    kernel_params = experiment_params['kernel']
    if dump_kernel_params:
        return neuron_params, kernel_params
    else:
        return neuron_params


def RunNetGrowth(n_samples, n_procs, neuron_params, save_path="tmp_measure",
                 plot=False):
    kernel = {"seeds": [33, 57, 19, 37, 79, 87, 11][:n_procs],
              "num_local_threads": n_procs,
              "environment_required": False}
    experiment_params = {}
    experiment_params["num_neurons"] = n_samples
    initial_position = [0, 0]
    positions = []
    for x in range(experiment_params['num_neurons']):
        positions.append(initial_position)
    neuron_params["position"] = positions
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())

    neuron_params['growth_cone_model'] = 'random_walk'

    NetGrowth.CreateNeurons(experiment_params["num_neurons"],
                                   "self_ref_forces",
                                   params=neuron_params,
                                   axon_params=neuron_params,
                                   dendrite_params=None,
                                   num_neurites=1,
                                   set_position=True)
    NetGrowth.Simulate(sim_length)
    if plot:
        NetGrowth.PlotNeuron()
    NetGrowth.SaveJson(filepath=save_path)
    NetGrowth.SaveSwc(filepath=save_path, swc_resolution=10)
    # NetGrowth.SaveJson(filepath=tmp_dir)
    NetGrowth.SaveSwc(filepath=os.path.join(
        os.getcwd(), save_path), swc_resolution=10)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()


def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)


def Test(neuron_params, plot=False):
    folder = os.path.join(os.getcwd(), "tmp_measure")
    CleanFolder(folder)
    RunNetGrowth(200, 5, neuron_params, folder)
    NG_population = NetGrowth.SimulationsFromFolder(folder)
    name = {"name": "self_ref_forces", "description":
            {"rw_memory_tau": 0,
             "rw_sensing_angle": 0,
             "rw_delta_corr": 0}
            }
    ensembles, fits = AnalyseNetgrowthRW(
        NG_population, int(max_len), info=name)
    # rw_corr.plot_results(ensembles, plot=True)
    info = InfoFromJson(os.path.join(folder, "info.json"))
    fit = fits[list(fits.keys())[0]]
    fit = analyse_fit(fit, info=info)
    fit = OnlyValues(fit)
    if plot:
        CleanFolder(folder)
        RunNetGrowth(1, 1, neuron_params, folder, True)
    CleanFolder(folder, make=False)
    print(" ################################### \n")
    # print(" Memory Tau: {} um \n".format(fit["memory"]))
    # print(" Correlation Tau: {} um \n".format(fit["pers_gauss"]))
    # print(" Sigma: {} \n".format(fit["sigma"]))
    print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    print(" ################################## \n")
    return fit["pers_length"], fit["tortuosity_local"], fit['cosine']

    # json.dump(fit,open(os.path.join(folder,"fit.json"),'w'))


if __name__ == "__main__":
    test = Test(neuron_params, True)

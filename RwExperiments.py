#!/usr/bin/env python
#-*- coding:utf-8 -*-

import NetGrowth
import numpy as np
import matplotlib.pyplot as plt
import os, sys

neuron_params = {
    "axon_angle":0.,
    "use_tubulin": False,
    "speed_growth_cone": 0.2,
    }

def step(n, loop_n, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        NetGrowth.PlotNeuron(show_nodes=True)

def automate(n_samples, n_procs, neuron_params, save_path, name):
    kernel={"seeds":[33, 57,19,37,79,87,11][:n_procs],
            "num_local_threads":n_procs,
            "environment_required":False}
    experiment_params={}
    experiment_params["num_neurons"]=n_samples
    initial_position=[400,-50]
    positions=[]
    for x in range(experiment_params['num_neurons']):
        positions.append(initial_position)
    neuron_params["position"]   =positions
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())

    neuron_params['growth_cone_model']='random_walk'

    gids =None
    plt.show()
    gids = NetGrowth.CreateNeurons(     n=experiment_params["num_neurons"],
                                        params=neuron_params,
                                        # dendrites_params=dendrite_params,
                                        num_neurites=1,
                                        set_position=True)
    NetGrowth.Simulate(10000)
    os.makedirs(os.getcwd()+"/"+save_path+"/"+name)
    NetGrowth.SaveJson(filepath=os.getcwd()+"/"+save_path+"/"+name)
    NetGrowth.SaveSwc (filepath=os.getcwd()+"/"+save_path+"/"+name,swc_resolution = 20)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()

##diffusion:
def reset_dict():
    neuron_params = {
        "axon_angle":0.,
        "use_tubulin": False,

    "rw_sensing_angle":0.1433,
    "speed_growth_cone": 1.,
        }
    return neuron_params

def test_diffusion(folder):
    exp_name=folder+"/diffusion_down/"
    neuron_params=reset_dict()
    os.makedirs(exp_name)
    for x in np.logspace(0,-2,4) :
        neuron_params["rw_sensing_angle"]=x
        neuron_params["rw_persistence_length"]= 0.000001
        automate(100, 5, neuron_params, exp_name, "sensing_angle_0"+str(x))

def test_memory(folder):
    neuron_params=reset_dict()
    neuron_params["rw_sensing_angle"]=0.1
    exp_name =folder+"/memory_without_pers_tris/"
    os.makedirs(exp_name)
    for x in np.logspace(1,3,10):
        neuron_params["rw_memory_tau"]=x
        neuron_params["rw_delta_corr"]=0.0001
        automate(100, 5, neuron_params, exp_name, "rw_memory_tau_0"+str(x))

def memory_vs_pers(folder,delta_corr):
    neuron_params=reset_dict()
    neuron_params["rw_sensing_angle"]=delta_corr
    exp_name=folder+"/memory_vs_pers__var_"+str(delta_corr)
    os.makedirs(exp_name)
    for y in np.logspace(1,3,10):
        for x in np.logspace(1,3,10):
            neuron_params["rw_memory_tau"]=x
            neuron_params["rw_delta_corr"]=y
            automate(100, 5, neuron_params, exp_name,"rw_memory_tau_0"+str(x)+"corr_tau"+str(y) )

def test_variance(folder):
    neuron_params=reset_dict()
    neuron_params["rw_persistence_length"]=1.
    exp_name=folder+"/test_variance/"
    os.makedirs(exp_name)
    for x in [ 100, 200, 500, 1000.]:
        neuron_params["rw_memory_tau"]=50.
        automate(x, 5, neuron_params, exp_name)

def linear_mem_pers(folder):
    neuron_params=reset_dict()
    exp_name=folder+"/linear/"
    os.makedirs(exp_name)
    for x in np.logspace(1,3,10):
        neuron_params["rw_persistence_length"]=x
        automate(200, 5, neuron_params, exp_name, "rw_pers_lenght_0"+str(x))

def quadratic_mem_pers(folder):
    neuron_params=reset_dict()
    exp_name=folder+"/linear/"
    os.makedirs(exp_name)
    for x in np.logspace(1,3,10):
        neuron_params["rw_delta_corr"]=x*x
        neuron_params["rw_memory_tau"]=x
        automate(200, 5, neuron_params, exp_name, "rw_memory_tau_0"+str(x))

def test_length(folder):
    neuron_params=reset_dict()
    exp_name=folder+"/linear/"
    os.makedirs(exp_name)
    neuron_params["rw_persistence_length"]=80.
    for x in np.logspace(-2,2,10):
        neuron_params["speed_growth_cone"]=x
        automate(200, 5, neuron_params, exp_name, "speed_0"+str(x))

def test_test(folder):
    neuron_params = {
        "axon_angle":0.,
        "use_tubulin": False,

    "rw_persistence_length": 0.,
    "rw_sensing_angle":0.1433,
    "speed_growth_cone": 1.,
        }
    exp_name=folder+"/test/"
    os.makedirs(exp_name)
    for x in [ 10., 20.]:
        neuron_params["rw_persistence_length"]=x
        automate(5, 5, neuron_params, exp_name)


folder = sys.argv[1]
test_diffusion(folder)
# test_memory(folder)
# test_variance(folder)
# linear_mem_pers(folder)
# memory_vs_pers(folder,0.5)
# memory_vs_pers(folder,1.0)
# memory_vs_pers(folder,100)

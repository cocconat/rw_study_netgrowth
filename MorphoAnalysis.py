#!/usr/bin/env python
#-*- coding:utf-8 -*-

# This program perform some measures on existing NetGrowth experiments
# it evaluates the random walk properties of a certain swc file.

# This software is part of NetGrowth project and SENEC initiative.

import os, json
from NetGrowthRwBenchmark import AnalyseNetgrowthRW, analyse_fit, plot_results, SwcToSegments, SegmentsToNetgrowth, fit_from_file, plot_fits, print_fits, OnlyValues
import NetGrowth
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate Random Walk properties for SWC files. It can be used to evaluate file generated through NetGrowth or downloaded by the NeuroMorpho.org archive.')
parser.add_argument('--fit', action='append', default =None,
                    help = 'Add `fit.json` file generated from the MorphAnlaysis to fit list to confront with other results, write --fit for each file')
parser.add_argument('--fit_parameter',type=str, help = "the parameter in respect to perform the fit")
parser.add_argument('--folder','-f', type=str, default=None,
                    help='Folder with NetGrowth experiment, if folder contains a `morphology.swc` file analyze it, other way analyze all the subfolders')
parser.add_argument('--neuron','-n', type=str, default=None, action= 'append',
                    help='Folder with SWC files downloaded from NeuroMorpho')
parser.add_argument('--max_len','-l', type=str, default=100,
                    help='max len, in micron, to analyze. 100 [um] is default. Micrometers is the standard unit for Swc files')

args = parser.parse_args()

def printhelp():
    print("the first argument is the folder with the simulation, \n"
            "the second argument is the max length to analyse")

if args.folder:
    """
    Analyze NetGrowth experiments
    """
    folder=args.folder
    if os.path.isdir(folder):
        print ("importing neurons from ", folder)
        NG_populations = NetGrowth.SimulationsFromFolder(folder)
    ensembles, fits =AnalyseNetgrowthRW(NG_populations,int(args.max_len))
    json.dump(fits,open(os.path.join(folder,"fit.json"),'w'))
    plot_results(ensembles, folder, plot=True)

if args.neuron:
    minimum_path=150
    max_angle=0.9
    names = ["axon", "dendrites"]
    Paths_populations=[]
    name=[ neuron.split("/")[1].split(".")[0] for neuron in args.neuron]
    Paths_populations.extend([np.concatenate([SwcToSegments(os.path.join(os.getcwd(),neuron), max_angle, minimum_path, element_type=[3]) for neuron in args.neuron])])
    NG_populations = [SegmentsToNetgrowth(paths, names[n], "experiment") for n,paths in enumerate(Paths_populations)]
    ensembles, fits =AnalyseNetgrowthRW(NG_populations,int(args.max_len))
    json.dump(fits,open("retina_fit.json",'w'))
    plot_results(ensembles, "retina_plot", plot=True)
    fit = OnlyValues(analyse_fit(fits['axon']))
    # dendrites = OnlyValues(fits['dendrites'])
    print(" ################################### \n")
    print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    print(" ################################## \n")

if args.fit is not None:
    fits_list=[]
    for fit_file in args.fit:
        fit = fit_from_file(fit_file)
        fits_list.append(analyse_fit(fit))
    print_fits(fits_list)
    plot_fits(fits_list, args.fit_parameter)


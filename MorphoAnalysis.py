#!/usr/bin/env python
#-*- coding:utf-8 -*-

# This program perform some measures on existing NetGrowth experiments
# it evaluates the random walk properties of a certain swc file.

# This software is part of NetGrowth project and SENEC initiative.

import os, json
from rw_corr import AnalyseNetgrowthRW, analyse_fit, plot_results, SwcToSegments, SegmentsToNetgrowth, fit_from_file, plot_fits, print_fits
import NetGrowth
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate Random Walk properties')
parser.add_argument('--fit', action='append',
                    help = 'add this file to fit list, as many as user reuired')
parser.add_argument('--fit_parameter',type=str, help = "which parameter in respect to perform the fit")
parser.add_argument('--folder','-f', type=str, default=None,
                    help='Folder with NetGrowth experiment')
parser.add_argument('--neuron','-n', type=str, default=None, action= 'append',
                    help='Folder with SWC file, from NeuroMorpho')
parser.add_argument('--max_len','-l', type=str, default=None,
                    help='max len to analyze. 100 [um] is fine.')

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

if args.fit !=[]:
    fits_list=[]
    for fit_file in args.fit:
        fit = fit_from_file(fit_file)
        fits_list.append(analyse_fit(fit))
    print_fits(fits_list)
    plot_fits(fits_list, args.fit_parameter)


import sys, os, json
from rw_corr import AnalyseNetgrowthRW, analyse_fit, plot_results, SwcToSegments, SegmentsToNetgrowth, fit_from_file, plot_fits, print_fits, plot_heatmap
import NetGrowth
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate Random Walk properties')
parser.add_argument('--mesure','-m', action='store_true', default=None,
                    help='experiment folder')
parser.add_argument('--fit', action='append')
parser.add_argument('--fit_parameter',type=str)
parser.add_argument('--folder','-f', type=str, default=None,
                    help='experiment folder')
parser.add_argument('--neuron','-n', type=str, default=None, action= 'append',
                    help='swc file from real neuron')
parser.add_argument('--max_len','-l', type=str, default=None,
                    help='maximum length to process')

args = parser.parse_args()

def printhelp():
    print("the first argument is the folder with the simulation, \n"
            "the second argument is the max length to analyse")

if len(sys.argv)<2 or sys.argv[1]=="--help":
    printhelp()
##
if args.folder:
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
    # Paths_populations = [np.concatenate([SwcToSegments(os.path.join(os.getcwd(),neuron), max_angle, minimum_path, element_type=[2]) for neuron in args.neuron])]
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


if args.heatmap:
    plot_heatmap(args.heatmap)


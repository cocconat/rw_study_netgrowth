from .ensemble import Ensemble
from .analysis import AnalyseNetgrowthRW
from .single_to_ensemble import  SwcToSegments, SegmentsToNetgrowth
from.fit_analisys import analyse_fit, plot_fits, fit_from_file, print_fits, OnlyValues
from .plot_results import plot_results

if __name__=="__main__":
    pass
    # exps = run_walkers(params)
    # #exps = variance(params)
    # plot_analysis(exps)

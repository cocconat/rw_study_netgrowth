from .ensemblerw import EnsembleRW
from .analysis import AnalyseNetgrowthRW
from.fit_analisys import analyse_fit, plot_fits, fit_from_file, print_fits, OnlyValues
from .plot_results import plot_results
import os, shutil

def CleanFolder(tmp_dir, make=True):
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if make:
        os.mkdir(tmp_dir)

if __name__=="__main__":
    pass
    # exps = run_walkers(params)
    # #exps = variance(params)
    # plot_analysis(exps)


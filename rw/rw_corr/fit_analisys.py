import json, os, collections, random
import sys
import numpy as np
import matplotlib.pyplot as plt

def fit_from_folder(fit_folder):
    jsons={json_file.split(".")[0]:json.load(open(fit_folder+json_file,'r')) for json_file in os.listdir(fit_folder) if json_file.split(".")[1]=="json"}
    return jsons

def fit_from_file(fit_file):
    jsons={fit_file.split(".")[0]:json.load(open(fit_file,'r'))}
    return jsons

def analyse_fit(jsons):
    """
    Merge all the file in the 'fit_folder' and plot results together,
    otherwhise plot the results for a certain 'fit_file'

    Params:
    fit_folder: folder with 'json' file from ensemble.fits
    fit_file:   'json' file from ensemble.fits

    Returns:
    results: dictionary with fit results.
    """
    def get_gauss(x):
        valuetext=x.split(" ")[0]
        value=valuetext.split("_")[1]
        return float(value)

    def get_variance(x):
        valuetext=x.split(" ")[2]
        value=valuetext.split("_")[1]
        return float(value)
    def get_memory(x):
        valuetext=x.split(" ")[1]
        value=valuetext.split("_")[1]
        return float(value)
    fits={}
    for key in jsons:
        for run in jsons[key]:
            fits[run]=jsons[key][run]
            ## get values from simulation:
            try:
                fits[run]['memory']={'values':{'a0':get_memory(run)},
                                     'errors':{'a0':0}}
                fits[run]['pers_gauss']={'values':{'a0': get_gauss(run)},
                                     'errors':{'a0':0}}
                fits[run]['variance'] = {'values':{'a0': get_variance(run)},
                                     'errors':{'a0':0}}
            except:
                fits[run]['name']=run

            ## get values from measures
            pers_error=fits[run]['msd_1D_ramp']['errors']['a0']/fits[run]['msd_1D_ramp']['values']['a0']
            fits[run]['pers_length']={'values':{'a0': 1./fits[run]['msd_1D_ramp']['values']['a0']},
                                     'errors':{'a0':pers_error/fits[run]['msd_1D_ramp']['values']['a0']}}
            fits[run]['tortuosity_dm']={'values':{'a0': 1./fits[run]['tortuosity_dm']['values']['a0']},
                                     'errors':{'a0':0}}

            if 'a2' in fits[run]['msd_2D']['values']:
                fits[run]['msd_2D']={'values':{'a0': fits[run]['msd_2D']['values']['a1']},
                                     'errors':{'a0':0}}
            fits[run]['name']=run
    for key in fits[run]:
        print (key)
    return fits

def memory_vs_pers(fits, x_a=0, y_a=0, ax=None, m='o', color=0) :
    x = np.linspace(1,3,10)
    y = np.linspace(1,3,10)
    xv, yv =np.meshgrid(x,y)
    z=[]
    for xx in x:
        for yy in y:
            for exp in fits:
                if fits[exp]['memory']['values']['a0'] == xx and fits[exp]['pers_gauss']['values']['a0'] == yy:
                    z.append(fits[exp]['pers_length']['values']['a0'])

    ax.pcolor(x,y,z)

    plt.loglog()
    # ax.legend()
    ax.set_xlabel("persistence")
    ax.set_ylabel("memory")
    return ax

def plot_experiment(fits, yvalue, x_a=0, y_a=0, ax=None, m='o', color=0) :
    if ax==None:
        fig,ax=plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None)
    c=[
    'b',
    'g',
    'r',
    'c',
    'm',
    'y',
    'k',
    'c',
    'm',
    'y',
    'k',
    'w'
    ]

    ax.scatter(range(len(fits.values())),
                  [fit[yvalue]['values']['a'+str(y_a)] for fit in fits.values()],
                  s=40,marker=m,color=c[color])
    ax.set_yscale('log')
    ax.set_xlabel("experiments")
    ax.set_ylabel(yvalue)
    return ax

def print_fits(fits_list):
    print(json.dumps(fits_list, sort_keys=True, indent=4))

def confront(fits, xvalue, yvalue, x_a=0, y_a=0, ax=None, m='o', color=0) :
    if ax==None:
        fig,ax=plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None)
    c=[
    'b',
    'g',
    'r',
    'c',
    'm',
    'y',
    'k',
    'c',
    'm',
    'y',
    'k',
    'w'
    ]

    ax.scatter([fit[xvalue]['values']['a'+str(x_a)] for fit in fits.values()],
                  [fit[yvalue]['values']['a'+str(y_a)] for fit in fits.values()],
                  s=40,marker=m,color=c[color])

    ax.errorbar([fit[xvalue]['values']['a'+str(x_a)] for fit in fits.values()],
                  [fit[yvalue]['values']['a'+str(y_a)] for fit in fits.values()],
                  [fit[yvalue]['errors']['a'+str(y_a)] for fit in fits.values()],
                  color=c[color])
    plt.loglog()
    ax.set_xlabel(xvalue)
    ax.set_ylabel(yvalue)
    return ax

# x= analyse_fit(fit_file="simulations/linear_fit.json")
# fits40=analyse_fit(fit_file="simulations/fits/memory_vs_pers40_fit.json")
# fits20=analyse_fit(fit_file="simulations/fits/memory_vs_pers20_fit.json")
# #fits10=analyse_fit(fit_file="random_walk_parameters/memory_vs_pers10_fit.json")
# #fits5  =analyse_fit(fit_file="random_walk_parameters/memory_vs_pers5_fit.json")
# fits    =analyse_fit(fit_file="simulations/fits/memory_fit.json")
# fits_list=[fits,fits20,fits40]
def plot_fits(fits_list, parameter):
    g, ((ax1, ax2), (ax3, ax4)  )  = plt.subplots(2,2,sharex=True)

    ax1.set_title("Persistence_length")
    ax3.set_title("Mean square displacement")
    ax2.set_title("Contraction ")
    ax4.set_title("Tortuosity local")

    # for n, fits in enumerate(fits_list):
        # for x in fits:
            # ax1=confront(x,"memory","pers_length",ax=ax1, color=n)
            # ax3=confront(x,"memory","msd_2D",ax=ax3, color=n)
            # ax2=confront(x,"memory","tortuosity_dm",ax=ax2, color=n)
            # ax4=confront(x,"memory","tortuosity_local",ax=ax4, color=n)

    for n, fits in enumerate(fits_list):
        if parameter=="name":
            ax1=plot_experiment(fits,"pers_length",ax=ax1, color=n)
            ax3=plot_experiment(fits,"msd_2D",ax=ax3, color=n)
            ax2=plot_experiment(fits,"tortuosity_dm",ax=ax2, color=n)
            ax4=plot_experiment(fits,"tortuosity_local",ax=ax4, color=n)
        else:
            ax1=confront(fits,parameter,"pers_length",ax=ax1, color=n)
            ax3=confront(fits,parameter,"msd_2D",ax=ax3, color=n)
            ax2=confront(fits,parameter,"tortuosity_dm",ax=ax2, color=n)
            ax4=confront(fits,parameter,"tortuosity_local",ax=ax4, color=n)
    ax2.legend()
    ax2.set_ylim([0.000001,0.001])
    ax2.semilogy( )
    ax1.semilogy( )
    ax3.semilogy( )
    ax4.semilogy( )
    g.tight_layout()
    g.show()
    g.savefig(os.getcwd()+"/fits.pdf",dpi=300,format="pdf")

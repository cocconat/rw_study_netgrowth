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
    # ax.errorbar([fit[xvalue]['values']['a'+str(x_a)]  for fit in fits.values()],
    #               [fit[yvalue]['values']['a'+str(y_a)] for fit in fits.values()],
    #               xerr=fit[xvalue]['errors']['a'+str(y_a)],
    #               yerr=fit[yvalue]['errors']['a'+str(y_a)], alpha=0.2, color=color)
    plt.loglog()
    # ax.legend()
    ax.set_xlabel(xvalue)
    ax.set_ylabel(yvalue)
    return ax

def plot_fits(fits_list, parameter):
    g, ((ax1, ax2), (ax3, ax4)  )  = plt.subplots(2,2,sharex=True)

    ax1.set_title("Persistence_length")
    ax3.set_title("Mean square displacement")
    ax2.set_title("Tortuosity R/L")
    ax4.set_title("Tortuosity local")

    # for n, fits in enumerate(fits_list):
        # for x in fits:
            # ax1=confront(x,"memory","pers_length",ax=ax1, color=n)
            # ax3=confront(x,"memory","msd_2D",ax=ax3, color=n)
            # ax2=confront(x,"memory","tortuosity_dm",ax=ax2, color=n)
            # ax4=confront(x,"memory","tortuosity_local",ax=ax4, color=n)

    for n, fits in enumerate(fits_list):
        for x in fits:
            ax1=confront_grid(x,parameter,"pers_length",ax=ax1, color=n)
            ax3=confront_grid(x,parameter,"msd_2D",ax=ax3, color=n)
            ax2=confront_grid(x,parameter,"tortuosity_dm",ax=ax2, color=n)
            ax4=confront_grid(x,parameter,"tortuosity_local",ax=ax4, color=n)
    ax2.legend()
    ax2.set_ylim([0.000001,0.001])
    ax2.semilogy( )
    ax1.semilogy( )
    ax3.semilogy( )
    ax4.semilogy( )
    g.tight_layout()
    g.savefig(os.getcwd()+"/fits.pdf",dpi=300,format="pdf")


def plot_heatmap(fit_file)

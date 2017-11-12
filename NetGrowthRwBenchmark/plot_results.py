import matplotlib.pyplot as plt
# from uncertainties import unumpy
import numpy as np
from matplotlib import colors as mcolors
import os

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names=[
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

def yerr(my_array):
    return [my_array[:,0]-my_array[:,2],my_array[:,1]-my_array[:,0]]

def plot_results(ensembles, save_path=None, plot=False):
    """
    Plot the path, the autocorrelation function, the variance, and the histogram
    for a set of 'experiments'.
    exps expect this format:
    [(path1, 'exp1_name'),...,(path2,'exp2_name')]
    """

    # f, axarr = plt.subplots(2)
    g, ((ax1, ax2), (ax3, ax4))  = plt.subplots(2,2,sharex=True)
    f,bx=plt.subplots()
    for n, ensemble in enumerate(ensembles):

        ax1.set_title("Contraction")
        ax1.set_xlabel("length um")
        ax1.set_ylabel(r"$\frac{Euclidean}{Length}$|")
        ax1.plot(ensemble.effective_length,
                ensemble.tortuosity_dm[:,0],
                label=ensemble.name,
                c=sorted_names[n])


        ax1.errorbar(ensemble.effective_length,
                    ensemble.tortuosity_dm[:,0],
                    c=sorted_names[n],
                    yerr=yerr(ensemble.tortuosity_dm),
                    alpha=0.5,
                    errorevery=10+n)

        ax1.plot(ensemble.effective_length,
                np.exp(- ensemble.effective_length / ensemble.fits['tortuosity_dm'][0][0])
                * ensemble.fits['tortuosity_dm'][0][1]+ensemble.fits['tortuosity_dm'][0][2],
                '--', c=sorted_names[n])

        ax2.set_title("position MSD")
        ax2.set_xlabel("length um")
        ax2.set_ylabel(r"$\langle X^2(l) \rangle$")
        ax2.plot(ensemble.effective_length,
                ensemble.msd_2D[:,0],
                label=ensemble.name,
                c=sorted_names[n]
                )
        ax2.errorbar(ensemble.effective_length,
                    ensemble.msd_2D[:,0],
                    yerr=yerr(ensemble.msd_2D),
                    c=sorted_names[n],
                    alpha=0.5,
                    errorevery=10+n)

        if ensemble.msd_2D_fit_wrong == 'msd_2D_quad':
        ##is linear, since the 2 fits are equal
            ax2.plot(ensemble.effective_length,
                    ensemble.effective_length
                    *ensemble.fits['msd_2D'][0][0] +
                    ensemble.fits['msd_2D'][0][1],
                    '--', c=sorted_names[n])
        else:
            ## is quadratic and the linear is wrong!
            ax2.plot(ensemble.effective_length  ,
                    ensemble.effective_length**2
                    *ensemble.fits['msd_2D'][0][0] +
                    ensemble.effective_length
                    *ensemble.fits['msd_2D'][0][1] +
                    # ensemble.effective_length*ensemble.fits['msd_2D'][0][1] +
                    ensemble.fits['msd_2D'][0][2], '--', c=sorted_names[n])

        ax3.set_title("Local tortuosity\n")
        ax3.set_xlabel("length um")
        ax3.set_ylabel("tortuosity")
        ax3.plot(ensemble.effective_length,
                ensemble.tortuosity_local[:,0],
                label=ensemble.name,
                c=sorted_names[n])

        ax3.errorbar(ensemble.effective_length,
                    ensemble.tortuosity_local[:,0],
                    yerr=yerr(ensemble.tortuosity_local),
                    c=sorted_names[n],
                    alpha=0.5,
                    errorevery=10+n)
        # ax3.plot(ensemble.effective_length,np.ones(len(ensemble.effective_length))*
        #     ensemble.fits['tortuosity_local'][0][0], '--', c=sorted_names[n], alpha=0.3)
        ax3.plot(ensemble.effective_length,\
                  np.ones(len(ensemble.effective_length))*ensemble.fits['tortuosity_local'][0][0], '--', c=sorted_names[n])

        # variance = [path[x]*path[x] for x in np.arange(10, len(path),10)]
        ax4.set_title("angle MSD")
        ax4.set_xlabel(" length um")
        ax4.set_ylabel(r"$\langle \theta^2(l)\rangle$")
        def get_transient(theta):
            for n,x in enumerate(theta):
                if x > 4:
                    break
            return n
        theta_max=get_transient(ensemble.msd_1D[:,0])
        ax4.plot(ensemble.effective_length,
                ensemble.msd_1D[:,0],
                label=ensemble.name,
                c=sorted_names[n])
        ax4.errorbar(ensemble.effective_length,
                ensemble.msd_1D[:,0],
                yerr=yerr(ensemble.msd_1D),
                c=sorted_names[n],
                alpha=0.3,
                errorevery=5+5*n)
        ax4.plot(ensemble.effective_length[0:theta_max],
                ensemble.effective_length[0:theta_max]*ensemble.fits['msd_1D_ramp'][0][0] +
                ensemble.fits['msd_1D_ramp'][0][1], '--',
                c=sorted_names[n])

    handles, labels = ax1.get_legend_handles_labels()
    bx.axis("off")
    bx.legend(handles, labels, loc='center')


    # ax1.legend(loc='right', shadow=True)
    g.tight_layout()
    # g.tight_layout()
    if save_path is not None:
        # f.savefig(os.getcwd()+"/"+save_path+"_positions_correlation.pdf",dpi=300,format="pdf")
        g.savefig(os.getcwd()+"/"+save_path+"_mean_square_disp__histogram.pdf",dpi=300,format="pdf")
        f.savefig(os.getcwd()+"/"+save_path+"_legend_.pdf",dpi=300,format="pdf")

    if plot:
        plt.show()

import json
from matplotlib import rc
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

z =[]

####
# This program compute <b_i b_j> for the memory constrained walker
# to do this it will compute <cos deltatheta> on an ensemble, and will measure for different \alpha.
# once compute the cosine it's possible to look at the correlation function and compute
# the persistence length


def correlate(path):
    """
    Return the correlation of path:
    C(i) = <path[i]*path[j]> / <path^2>
    """
    path=np.array(path)
    z=[]
    var = np.mean(path*path)
    for n in range(1,path.shape[0]):
        z.append(np.mean(path[:-n]*path[n:])/var)
    return np.array(z)

def cosine_delta(array, n):
    """
    Return cosine of delta
    """
    z = np.array(array[n:]) - np.array(array[n])
    return np.cos(z)


def Run(time, alpha, sigma):
    """
    Run the process with memory factor 'alpha' and 'time' timestep

    Return
    ======
    xy: tuple of np.array
    theta: np.array
    """
    theta = 0
    theta_bar = 0
    pathx =[0]
    pathy =[0]
    path =[]
    alpha_temp_sum =0
    for x in range(time):
        theta_bar = (theta + alpha_temp_sum * alpha * theta_bar )/ (alpha * alpha_temp_sum +1)
        alpha_temp_sum = alpha *alpha_temp_sum +1
        theta = theta_bar + (np.random.rand()-0.5)*sigma
        pathx.append(pathx[-1]+np.cos(theta))
        pathy.append(pathy[-1]+np.sin(theta))
        path.append(theta)

    return (pathx,pathy), path

def SmoothSlide(path, l=10):
    """
    Smooth sequence with slide window
    """
    out = np.copy(path)
    for n in range(out.shape[0]):
        out[n]=np.mean(path[n:n+l])
    return out


def SquareDistance(xy):
    xy = np.array(xy)
    Xsquare = xy[0][1:]**2 + xy[1][1:]**2
    return Xsquare


def DeltaAutocorrelation(alpha, ensemble_size, time, sigma):
    """
    Compute the autocorrelation on Delta Theta angle.
    <Delta Theta (n)> over an ensemble of size 'ensemble_size'
    and with experiment parameters 'alpha' and 'time'

    Return
    ======
    delta = np.array(time)
    X^2 / X_0 ^2 where X_0 is diffusive process
    """
    z = np.ndarray((ensemble_size,time))
    k = np.ndarray((ensemble_size,time))
    for n in range(ensemble_size):
        xy, theta = Run(time, alpha, sigma)
        z[n] = cosine_delta(theta,0)
        k[n] = SquareDistance(xy)
    return np.mean(z, axis = 0) , np.mean(k, axis=0)/np.arange(1,time+1)

def FitExpDiverge(xdata, ydata, _plt=None):
    """
    Fit data with exponential diverging function

    Return
    ======
    m, q

    """
    def exp(x, l, q, x_0):
        return (q * np.exp(x *x * l))
    popt, pcov = optimize.curve_fit(exp, xdata, ydata, bounds =(0,100))
    # if _plt is None:
        # _plt = plt
    # _plt.plot(xdata, exp(xdata,*popt))
    # _plt.plot(xdata, ydata)
    return popt

def FitExpDecay(xdata, ydata, _plt=None):
    """
    Fit data with exponential decay

    Return
    ======
    m, q

    """
    def exp_low(x,l,q):
        return (q * np.exp(-x / l))
    popt, pcov = optimize.curve_fit(exp_low, xdata, ydata)
    # if _plt is None:
        # _plt = plt
    # _plt.plot(xdata, exp_low(xdata,*popt))
    # _plt.plot(xdata, ydata)
    return popt[0], popt[1]

def MeasureTau(data):
    data = SmoothSlide(data)
    return np.where(data<np.exp(-1))[0][0]



def run_alpha(sigma, alphas, _time, ensemble_size):
    """
    Run the experiment for several alpha coefficient and return a dictionary with all results.
    """
    out ={}
    for alpha in alphas:
        time = int(_time *alpha)
        _delta, r_r0 =DeltaAutocorrelation(alpha, ensemble_size, time, sigma)
        _delta = SmoothSlide(_delta, 5)
        l_p, x_0 = FitExpDecay(np.arange(time), _delta )
        ### Here we compute the volume below the distribution S,
        ### and the first moment M
        ### for exponential function a*exp(-l/l_p) the integral S is a*l_p
        ### the first moment M is l_p**2
        S=x_0*l_p
        M=l_p**2*x_0
        out[alpha] = {"l_p":l_p, "x_0":x_0, "S":S, "M":M, "R/R0": r_r0[-1]}
        ax1.plot(_delta,label = alpha)
        ax2.scatter(alpha, l_p)
        ax4.plot(r_r0/r_r0[-1])
        ax3.scatter(alpha, (1 + 2*(S -M/time)) , c='r',m='*', label="S")
        ax3.scatter(alpha, r_r0[-1] , c='b',m='.', label="R_R0")
    lps = [out[x]["l_p"] for x in alphas]
    tau_alpha = 0
    q_0_alpha = 0
    # tau_alpha, q_0_alpha, x_0= FitExpDiverge(alphas, lps, ax2)
    return { "out":out, "tau_alpha": tau_alpha, "q_0": q_0_alpha, "x_0":x_0}



#####################################################
### Run the experiment.
#####################################################
ensemble_size =  1000
_time = 10000

fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(2,2)

alphas = np.linspace(0.2,0.9,8)
sigma_dict={}
sigmas=[0.05, 0.1, 0.5, 1., 1.5]
sigmas=[1]

for sigma in sigmas :
    sigma_dict[sigma] = run_alpha(sigma, alphas, _time, ensemble_size)

# ax1.scatter(list(sigma_dict.keys()), [sigma_dict[key]["tau_alpha"] for key in sigma_dict.keys()])
# ax2.scatter(list(sigma_dict.keys()), [sigma_dict[key]["q_0"] for key in sigma_dict.keys()])
# ax3.scatter(list(sigma_dict.keys()), [sigma_dict[key]["x_0"] for key in sigma_dict.keys()])


ax1.legend()
ax1.set_title(r'$ \langle b_i b_j \rangle $')
ax2.set_title(r'Persistence Length')
ax2.set_title(r'Persistence Length')
ax1.set_xlabel(r'time')
ax1.set_ylabel(r'$\cos(\Delta \theta)$' )
ax2.set_ylabel(r'$l_p$')
ax2.set_xlabel(r'$\alpha$')
ax3.set_title(r'S and R/R0')
ax4.set_title(r'Convergence to diffusive behavior')
ax4.set_xlabel(r'time')
ax3.set_xlabel(r'$\alpha$')
ax4.set_ylabel(r'R/R0')

plt.tight_layout()
plt.savefig("Volume_memory_study3.pdf", dpi = 300)

# with open ("Volume_memory",'w') as fp:
    # json.dump(out,fp)



plt.show()
# plt.plot(np.exp(-np.arange(0,1000)/1000)-0.3)


    # memory_.effective_angle =
            # (move_.angle +
             # memory_.alpha_temp_sum * memory_.alpha_coeff *
                 # memory_.effective_angle) /
            # (memory_.alpha_coeff * memory_.alpha_temp_sum + 1);

        # memory_.alpha_temp_sum =
            # memory_.alpha_coeff * memory_.alpha_temp_sum + 1;
    # move_.angle = memory_.effective_angle;


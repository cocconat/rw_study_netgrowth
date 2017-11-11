from .random_walks import correlated_rw, classic_rw, correlated_gaussian, initialize_exp
from .algorithms import msd_1D, msd_2D, tortuosity_local, tortuosity_dm
from NetGrowth import get_axon_path, get_properties
import numpy as np
from scipy import optimize
from scipy import signal as sgl
from uncertainties import unumpy

class Ensemble(object):

    def __init__(self, description):
        try:
            self.name, self.description =get_properties(description)
        except:
            self.name, self.description = description["name"], description["description"]
        self.theta= []
        self.r    = []
        self.xy   = []

    def add_sequence(self,neurons):
        for neuron in neurons:
            xy, polar = get_axon_path(neuron=neurons[neuron]['data'])
            self.theta.append(polar[1])
            self.r.append(polar[0])
            self.xy.append(xy)
        self.theta = np.array(self.theta)
        self.r     = np.array(self.r)
        self.xy    = np.array(self.xy)


    def characterize(self, max_len, first=1):
        """
        Run characterization algorithm on the random_walk ensemble:
        tortuosity_dm  : measure the tortuosity length over distance ratio between curvilinear distance and euclidean
        tortuosity_local : measure the correlation between successive delta.
        msd          : measure the mean square displacement of angles
        space_msd    : measure the mean square displacement of cartesian coordinates.

        all the measure are performed with a single istance of random walk, the average is perfermed
        over block with a distance longer than the correlation length.

        Params:
        ------
        max_len: maximal relative distance to measure for the different algorithms.
        decorrelation_length: distance of correlation to resize the path in uncorrelated segments.
        """
        if first < 1:
            raise Exception("first value has to be greater than one or dividion by zero occurs")
        self.max_len=max_len
        self.r=np.array(self.r)
        self.theta=np.array(self.theta)
        self.xy=np.array(self.xy)

        print ("max length is", self.max_len)
        self.tortuosity_local, self.max_len =  tortuosity_local(self.theta,self.r,self.max_len, first)
        print ("max length is", self.max_len)
        self.msd_1D= msd_1D(self.theta,self.max_len,first)
        self.msd_2D  = msd_2D(self.xy,self.max_len,first)
        self.tortuosity_dm  = tortuosity_dm(self.r,self.xy, self.max_len,first)
        self.effective_length=np.zeros((self.max_len))
        length=np.mean(self.r,axis=0)
        for shift in range(first,self.max_len+first):
            self.effective_length[shift-first]=np.sum(length[first:shift])

    def fit(self,last,first=0):
        def constant(x,k):
            return k
        def linear(x,m,b):
            return x*m+b
        def quadratic(x,a,b,c):
            return x**2*a+x*b+c
        def exponential(x, tau, a, b):
            return a*np.exp(-x/tau)+b
        def get_tau(array):
            array - 0.367879441
            for n,x in enumerate(array):
                if x<0:
                    return n
        def get_transient(theta):
            for n,x in enumerate(theta):
                if x > 2:
                    break
            return n

        self.fits={}
        # print (self.msd_1D[:,0])
        theta_max = get_transient(self.msd_1D[:,0])
        self.fits['msd_1D_ramp'] = optimize.curve_fit(linear,self.effective_length[:theta_max],
                self.msd_1D[:theta_max,0], sigma=np.abs(self.msd_1D[:theta_max,1]-self.msd_1D[:theta_max,2]))
        self.fits['msd_2D_quad'] = optimize.curve_fit(quadratic,self.effective_length,
                                        self.msd_2D[:,0])
        self.fits['msd_2D_lin'] = optimize.curve_fit(linear,self.effective_length,
                                        self.msd_2D[:,0])
        self.fits['tortuosity_dm'] = optimize.curve_fit(exponential,self.effective_length,
                self.tortuosity_dm[:,0], check_finite=True, sigma = np.abs(self.tortuosity_dm[:,1]-self.tortuosity_dm[:,2]), bounds=([100,0,0],[np.inf,1,1]))
        self.fits['tortuosity_local'] = optimize.curve_fit(constant,self.effective_length,
                                        self.tortuosity_local[:,0])
        if np.sum(np.diag(self.fits['msd_2D_quad'][1])**2) > np.sum(np.diag(self.fits['msd_2D_lin'][1])**2):
            self.msd_2D_fit_wrong='msd_2D_quad'
            self.fits.pop(self.msd_2D_fit_wrong)
            self.fits['msd_2D']= self.fits['msd_2D_lin']
            self.fits.pop('msd_2D_lin')
        else:
            self.msd_2D_fit_wrong='msd_2D_lin'
            self.fits.pop(self.msd_2D_fit_wrong)
            self.fits['msd_2D']= self.fits['msd_2D_quad']
            self.fits.pop('msd_2D_quad')
        self.results={}
        for key in self.fits:
            self.results[key]={"values":{"a"+str(n) : x for n, x in enumerate (self.fits[key][0])},
                             "errors":{ "a"+str(n) :x for n, x in enumerate (np.diag(self.fits[key][1]))}}
        return self.results

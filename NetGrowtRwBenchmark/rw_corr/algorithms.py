import numpy as np
from scipy import signal
import uncertainties as un
from uncertainties import unumpy



#============================
# Correlation functions
#============================


def tortuosity_local(theta,rho,max_len,first=1):
    """
    Measure the local tortuosity as defined in 'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    Params:
    -------
    arrays: numpy 2D array, first dimension is the set of realizations
    max_len: max interval where correlation is measured.
    x_0   : first element of the list.

    Returns:
    --------
    autocorrelation: 2D array, the percentile distribution with [50%, 75%, 25%]
    """

    # theta =  (arrays[:,:].transpose()- arrays[:,0]).transpose()
    # differential measure reduce size -1
    print(theta.shape)
    theta = np.abs( theta[:,1:] - theta[:,:-1])
    rho   =rho[:,1:]
    if max_len> len(rho[1,:-1]):
        print ("correcting max length to: ", len(rho[1,:-1])-10)
        max_len=len(rho[1,:-1])-10
    tortuosity_local = np.zeros((max_len,3))
    length_local = np.zeros((max_len,3))
    #since distance needs to be greater thean 1, first element is jumped!
    for shift in range(first+1,max_len+first):
        tortuosity_local[shift-first] = np.percentile(np.sum(theta[:,first:shift], axis=1),
                                                      q=[50,75,25], axis=0, interpolation='midpoint')
        # print(first, shift+1, rho[:,first:shift+1])
        length_local[shift-first] = np.percentile(np.sum(rho[:,first:shift], axis=1),q=[50,75,25], axis=0, interpolation='midpoint')
    length_local[0]=1,1,1
    tortuosity_local = tortuosity_local/length_local
    return tortuosity_local,len(tortuosity_local)

def tortuosity_dm(curvilinear, xy, max_len, first=1):
    """
    Measure the tortuosity_dm ratio between two points
    The tortuosity_dm is the ratio between the euclidean distance and curvilinear distance
    This value is always less than 1.

       Params:
    -------
    array: numpy 1D array
    max_len: double max interval where correlation is measured.

    Returns:
    --------
    autocorrelation: 2D array, the percentile distribution with [50%, 75%, 25%]
    """
    x=(xy[:,0,:].transpose() - xy[:,0,0]).transpose()
    y=(xy[:,1,:].transpose() - xy[:,1,0]).transpose()
    ratio =np.zeros((max_len,3))
    for shift in range(first,max_len+first):
        r  =np.sqrt((x[:,shift])**2 + (y[:,shift])**2)
        length = np.sum(np.abs(curvilinear[:,:shift]), axis=1)
        ratio[shift-first] = np.percentile(r/length, q=[50,75,25], interpolation='midpoint')
    return ratio

def msd_1D(array, max_len, first=1):
    """
    Compute the mean square displacement with a temporal average over non stationary sequence.
    Using the translational invariance of msd on such sequence
    Delta_x(n) = <(x_n - x_0)^2 > would require a set of replica to be averaged, we can compute
    Delta_x(m-n) = <(x_n - x_m)^2 > just ensuring the two different blocks are uncorrelated.
    The decorrelation time needs to be set and depends from the system

       Params:
    -------
    array: numpy 1D array
    max_len: double max interval where correlation is measured.
    decorrelate: time lapse to consider the replica indipendent each other

    Returns:
    --------
    msd: 2D array, the percentile distribution with [50%, 75%, 25%]
    """

    theta =  (array[:,:].transpose()- array[:,first]).transpose()
    msd = np.zeros((max_len,3))
    for shift in range(first,max_len+first):
        msd[shift-first] = np.percentile(theta[:,shift]**2,q=50, axis=0, interpolation='midpoint'),\
            np.percentile(theta[:,shift]**2,q=75, axis=0, interpolation='midpoint'),\
            np.percentile(theta[:,shift]**2, q=25, axis=0, interpolation='midpoint')
    msd[0][1] = 1
    msd[0][2] = 0
    return msd

def msd_2D(xy, max_len, first=1):
    """
    Compute the mean square displacement with a temporal average over non stationary sequence.
    Using the translational invariance of msd on such sequence
    Delta_x(n) = <(x_n - x_0)^2 > would require a set of replica to be averaged, we can compute
    Delta_x(m-n) = <(x_n - x_m)^2 > just ensuring the two different blocks are uncorrelated.
    The decorrelation time needs to be set and depends from the system

       Params:
    -------
    xy: numpy 2D array
    max_len: double max interval where correlation is measured.
    decorrelate: time lapse to consider the replica indipendent each other

    Returns:
    --------
    msd: numpy 2D array, the percentile distribution with [50%, 75%, 25%]
    """

    x=(xy[:,0,:].transpose() - xy[:,0,first]).transpose()
    y=(xy[:,1,:].transpose() - xy[:,1,first]).transpose()
    msd = np.zeros((max_len,3))
    for shift in range(first,max_len+first):
        msd[shift-first] = np.percentile(x[:,shift]**2+y[:,shift]**2, q=[50,75,25], axis=0, interpolation='midpoint')
    return msd

def msd_fft(r):
    """
    Compute the mean square displacement for the sequence
    """
    def autocorrFFT(x):
        N=len(x)
        F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
        PSD = F * F.conjugate()
        res = np.fft.ifft(PSD)
        res= (res[:N]).real   #now we have the autocorrelation in convention B
        n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
        return res/n #this is the autocorrelation in convention A
    N=len(r)
    if len(r.shape)>1:
        D=np.square(r).sum(axis=1)
        D=np.append(D,0)
        S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    else:
        D=np.square(r)
        D=np.append(D,0)
        S2=autocorrFFT(r)
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    return S1-2*S2

import btmorph2
import numpy as np
import os
from os.path import join, isfile
from os import listdir
import json

def get_properties(neuron):
        props = neuron['neurons']['0']['dendrites_params']
        name= "lp_{} mem_{} var_{}".format(props['rw_persistence_length'],
                                    props['rw_memory_tau'],
                                    props['rw_sensing_angle'])
        return name, props

def swc_to_array(neuron):
    return np.loadtxt(neuron, usecols=(2,3))

def import_swc(swc_path, population_to_singles=False):
    """
    Import the Neuron Morphology from swc.
    If the file contains more than a neuron and population_to_singles set to False
    creates a PopulationMorphology object,
    If the population_to_singles is True it always return a list of NeuronMorphology.

    Returns
    -------
    a list of bitmorph objects.
    """
    swc_file=swc_path+".swc"
    gids = btmorph2.neurons_from_swc(swc_file)
    neurons=[]
    for gid in range(1,gids+1):
        print("read {}/neuron_{}.swc".format(swc_path,gid))
        neurons.append(swc_to_array(swc_path+"/neuron_"+str(gid)+".swc"))
    return neurons, gids

def neurons_from_folder(folder_path, population_to_singles=False):
    """
    Return a list of netgrowth_format neurons for all the neurons in the folder.
    The folder is expected to be a set of .swc and .json files with same name.
    The .json file will contain information on the neurons.
    In case the swc file contain more than a neuron, let'say N,
    it will be splitted in N .swc files inside a folder with same name and path.
    This is done for compatibility with btmorph2.

    Returns
    -------
    [...,btmorph2.NeuronMorphology objects,...]
    """
    neuron_folder =os.getcwd()+"/"+folder_path+"/"
    print ("neuron_folder: ",neuron_folder)
    neuronfiles = [join(neuron_folder,f.split(".")[0]) for f in listdir(neuron_folder) if isfile(join(neuron_folder, f)) and f.endswith("swc")]
    neurons=[]
    for neuron in neuronfiles:
        print( "importing population from {}".format(neuron))
        imported_file, gids = import_swc(neuron, population_to_singles)
        print( "This population has {} neurons".format(gids))
        # print("imported_file: ", imported_file)
        netgrowth_format = {"gids":gids,"swc":imported_file,"json":json.load(open(neuron+".json"))}
        neurons.append(netgrowth_format)
    # neuronfiles = tuple_from_files(neuronfiles)
    return neurons

def get_neuron_path(neuron, plot=False):
    """
    Magic function!
    it recognizes the input format for neuron:
    * file path to swc file
    * btmorph NeuronMorphology object
    * np.ndarray
    and converts xy lists to xy and polar coordinates
    """
    if isinstance(neuron,np.ndarray):
        xy= neuron.transpose()
    elif isfile(neuron):
        neuron = btmorph2.neuron_from_file(neuron)
        neuron = neuron['swc']
    if isinstance(neuron, btmorph2.NeuronMorphology):
        if plot:
            neuron.plot_1D()
        xy=btmorph2.get_neuron_path(neuron)[:,5:]
    angles=angles_from_xy(xy)
    modules=module_from_xy(xy)
    return xy, np.array([modules,angles])

def angles_from_xy(path):
    angles=[]
    for n in range(1,len(path[0])):
        deltax=path[0][n]-path[0][n-1]
        deltay=path[1][n]-path[1][n-1]
        rad = np.arctan2(deltay,deltax)
        angles.append(rad)
    return demodularize(np.array(angles)-angles[0])

def demodularize(angles):
    shift=0
    demodule=np.zeros((len(angles)))
    for n, theta in enumerate(angles):
        if abs(theta-angles[n-1]) > 3.14:
            shift+=np.sign(angles[n-1])*2*np.pi
        demodule[n]=theta+shift
    return demodule

# def remove_modulus(angle, previous):
    # if angle
def module_from_xy(path):
    modules=[]
    for n in range(1,len(path[0])):
        deltax=path[0][n]-path[0][n-1]
        deltay=path[1][n]-path[1][n-1]
        module = np.sqrt(deltay**2 + deltax**2)
        modules.append(module)
    return modules

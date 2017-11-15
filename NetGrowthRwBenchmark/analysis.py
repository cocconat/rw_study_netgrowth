#!/usr/bin/env python
#-*- coding:utf-8 -*-
# This software is part of the NetGrowth project and the SENEC initiative

from . import Ensemble

def AnalyseNetgrowthRW(NG_populations, max_len=880, first=10):
    """
    Apply measures to all the neurons.swc in the folder, not the subfolder.
    Measures are those defined inside 'Ensemble.characterize',
    The measures are generally statistical measures computed over ensemble average,
    The ensemble is the set of neurons present in each .swc file in 'folder'.

    Params:
    ------

    folder : (str required) the folder with the 'swc','json' couple as saved from NetGrowth
    last: (int) the index of last element to average over the ensemble
    first: (int) the index of first element to average over the ensemble

    Returns:
    -------
    ensembles: list of Ensemble object
    fits: results from characterization

    """
    ensembles=[]
    fits={}
    for population in NG_populations:
        ensemble=Ensemble(population['info'])
        ensemble.add_population(population['neurons'])
        ensemble.characterize(max_len,first)
        fits[ensemble.name]=ensemble.fit(max_len,first)
        ensembles.append(ensemble)
    return ensembles, fits


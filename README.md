NetGrowthRWBenchmark
=====================


1. Installation
Every file is compiled on fly, then after downloading the repository from git to run the program you can just run the specific scripts

2. Scripts
This module offers 4 frontend scripts and a bunch of backend functions.
The frontend files are:
    - MeasurePersistence.py it's measure the `persistence length` and the `tortuosity` of the Negrowth random walker parameters set by the user. The parameters can be set at the beginning of the file. Persistence is measured in 2 manners, they are accurately: with MSD slope and cosine fitting of <b_ib_j>.
    - MemoryRwSimulator.py compute an ensemble of Memory Random Walk and compute it's properties (see above) and perform a fit. The file is quite self-explicative when you read it,
    - MorphoAnalysis.py is the tool for analyze already generated SWC files. Calling the command with `--help` option offers the whole set of operation you can do.
    - GridSearch.py searches for the values of `memory`, `sigma` and `delta_corr` required to obtain a plausible persistence length. It's required to set a range, with max and min persistence length at the beginning of the file.
    - RwExperiments.py create a set of NetGrowth experiment with the inserted properties, it's obsolete, use MeasurePersistence.py now.

3. BackEnd.
The idea is:
    1. Get a 'morphology.swc' file with N neurons inside and decompose it in N files, NetGrowth does it into dataIO_rw.py
    2. Import this files into `Ensemble` object. `analysis.py` and `ensemble.py`
    3. Ensemble characterizes the popultaion with tha algorithms in `algorithms.py`
    4. The fit are processed in `fit_analysis.py`

`single_to_ensemble.py` is required to apply statistical measures over a single neuron, this is done decomposing each neuron in many segments, cut at the branching points.

`rw_correlated.py` and `random_walk.py` are old proof of concept, ow correctly implemented into NetGrowth.

`plot_results.py` let plot the measures into 4 graphs scheme, tortuosity isn't plotted since it's a constant value.






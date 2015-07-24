#Python and C tools for studying Dynamical Systems

##Requirements
These tools require a C compiler, Python (version >= 2.7), and the following Python modules:
- [numpy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [graph_tool](https://graph-tool.skewed.de/)

##Installation and Building
I recommend installing the first two Python modules with the Python Package Index [pip](https://pypi.python.org/pypi/pip).
For graph\_tool, see [here](https://graph-tool.skewed.de/download). For OS X users we recommend using homebrew to install graph\_tool. We also caution that graph\_tool is large, and sometimes quite difficult to install. As it is only used in our Morse Decomposition program, the user may find graph\_tool unnecessary.

If you don't have a C compiler, we again recommend using homebrew (OS X) or apt-get (Linux) to get one. To build, simply type:

```
cd chaos
make -f make_all
```

##Tools
For the sake of speed, each python tool "xxrun.py" uses one or more compiled C libraries through Python's ctypes module. Included are the following tools:
 - A C Program to run 4th Order Runge Kutta in 2 Dimensions (rk4/rk4.c)
 - A C implementation of a Sparse Matrix for efficient graph storage (sparse\_matrix\_table/smtable.c)
 - A C Program with a corrected real number modulus (usefulfunctions/functions.c)
 - A C Program to compute the probability of Cantori switches (invariantdynamics/switch.c)
 - A Python Program to compute Invariant Sets (invariantdynamics/irun.py)
 - A Python Program to compute Birkhoff Partitions (invariantdynamics/bprun.py)
 - A Python Program to differentiate Quasiperiodic and Chaotic Sets (invariantdynamics/qrun.py)
 - A Python Program to compute Recurrence Plots (invariantdynamics/rprun.py)
 - A Python Program to compute Stickiness Rates along trajectories (invariantdynamics/srun.py)
 - A Python Program to compute Morse Decompositions of Invariant Sets (invariantdynamics/morse_analysis.py)
 - A Python Program to compute Sticky Regions in Chaotic Sets (invariantdynamics/sticky_analysis.py)
 - A Python Program to track the evolution of a region (figurefactory/frun.py)

##Usage
The parameters of each python program are explained in their comments. If the user wishes to add custom maps or flows, they must alter the corresponding C programs. Namely,
 - irun.py calls invariant4.c and smtable.c
 - bprun.py calls birkhoff_partition.c
 - qrun.py calls quasiperiodicity.c
 - rprun.py calls stick.c
 - srun.py calls sticky.c
 - rk4.c contains the differential equation governing any flow
To use one of the python programs, simply alter it's parameters as needs be and type
```
python xxrun.py
```
The two analysis programs don't have easily editable parameters, 
but are well commented if the user wishes to modify them. 

irun.py can also produce a graph file ('outputs/graph_save.p') if the user runs
```
python irun.py -g
```
This output is needed by morse\_analysis.py (don't run morse\_analysis.py unless you have produced a graph).
Similarly, sticky\_analysis.py requires the outputs of qrun.py and srun.py (run these first).

Finally, the program frun.py requires a directory named for the map being studied, 
and for the map to be specified in factory.c.

##Examples Website

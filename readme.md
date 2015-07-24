#Python and C tools for studying Dynamical Systems

##Requirements
These tools require a C compiler, Python (version >= 2.7), and the following python modules:
- numpy
- matplotlib
- ctypes
- cPickle
- graph_tool

##Installation and Building
I recommend installing the first four python modules with the Python Package Index (pip)[https://pypi.python.org/pypi/pip].
For graph\_tool, see (here)[https://graph-tool.skewed.de/download]. For OS X users we recommend using homebrew to install graph\_tool. We also caution that graph\_tool is large, and sometimes quite difficult to install. As it is only used in our Morse Decomposition program, the user may find graph\_tool unnecessary.

If you don't have a c compiler, we again recommend using homebrew (OS X) or apt-get (Linux). To build, simply type:

'''
cd chaos
make -f make_all
'''



##Tools

##Usage

##Website

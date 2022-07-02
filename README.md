# tomplot
The Only Matplotlib PDE Library Of Tom

tomplot is an interface to matplotlib for quickly making neat plots for
numerical models of PDEs.

It is intended for use with LFRic and some Firedrake packages, such as
transportdrake.

### Why?

tomplot is motivated by the importance of neat data visualisation. To avoid code
duplication across different projects or even plotting scripts, tomplot contains
a set of routines that are commonly used to make figures ready for publication
or presentation. This is all provided as an interface to matplotlib, so that all
these plots can be made quickly and easily, while allowing the full capability
of matplotlib to be exploited.

tomplot provides:
- routines to extract and process data
- default formats for various types of plot, to quickly make neat plots
- several specialised routines used in making publication-worthy plots

### INSTALLATION

The main dependencies are netCDF4, numpy and matplotlib. To take advantage of
full features, it requires pandas, scipy and cartopy.

Once downloaded, can be installed in a virtual environment with
```
python3 setup.py install
```
or made available through
```
export PYTHONPATH=<path_to_tomplot>/tomplot
```
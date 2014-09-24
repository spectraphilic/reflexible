# pflexible #

## About ##
This module is being created to work with FLEXPART model ouput. FLEXPART is a 
Lagrangian particle dispersion model written by Andreas Stohl. It is in 
active development at NILU and by Sabine Eckhardt, John Burkhart and numerous 
other research groups worldwide.

This module is based partially on Matlab(c) routines developed
by Sabine Eckhardt and Harald Sodemann at NILU. Subsequent development
and the Python routines have been conducted by John Burkhart.

This software was written for the purposes of conducting research and displaying
model results for publications. More information can be found on my
[homepage](http://niflheim.nilu.no/~burkhart) under the 'Software' and
'Publications' sections.

## Contact, Issues, and Discussion ##

A maillist has been set up for users at sourceforge. Please direct inquiries via
the [pflexible_maillist](https://lists.sourceforge.net/lists/listinfo/pflexible-users)

As for Issues, please use the
[tool](https://bitbucket.org/jfburkhart/pflexible/issues/new) at
bitbucket.

Direct inquiries may be sent to::

	John F Burkhart, University of Oslo
	[john.burkhart](mailto://john.burkhart@geo.uio.no)


## Web Resources ##

### BitBucket ###

Clone the bitbucket repository into a folder in your PYTHONPATH::

	%hg clone https://bitbucket.org/jfburkhart/pflexible pflexible

### Sphinx Documentation ###

The documentation is hosted at my
[homepage](http://niflheim.nilu.no/jfburkhart). You can get to it under the
'Software' section.

I am also trying to build it at readthedocs, but there are currently problems
with the autodoc for the modules [readthedocs](http://pflexible.readthedocs.org)


## Working with pflexible

There are a few 'gotchas' when using the module. First, you will likely
have to recompile (f2py) the FortFlex.f file and create a FortFlex.so module
for whatever computer you're using. The one presently is compiled for 
64bit Linux (Ubuntu 12.04). For more information see the f2py directory.

An alternative 'BinaryFile' class has been created so one can work with pure
Python. Alone, it is significantly slower than the FortFlex module, however,
if you use the dumpgrid module, significant speedups can be achieved. For this
run::

  python setup.py build_ext --inplace


This will compile the pflexcy.pyx file into a pflexcy.so module that can be
imported and used by the pf.readgrid function. A series of tests are run to try
and determine the best module to use -- somewhat transparently to the user. See
the pflexible.readgridV8 function for more information.

Primary functionality comes from the readheaderV8 function and the readgridV8
function. I have tried to create a "Header" class that can be used for some 
typical analysis. See the examples directory, and don't forget to read the
source code.

## Key Tools

The primary workflow and usefulness of pflexible comes from the read_header and read_grid
routines, which are designed to help ultimately with plotting and data analysis of
the FLEXPART output.

To get started, at least become familiar with:

* Header class
* read_header
* read_grid
* fill_backward
* plot_sensitivity

NOTE: the other plotting functions are mostly wrappers which
ultimately make a call to plot_sensitivity.

## Disclaimer

About the mapping.py module, it is primarily designed as some convenience 
routines for working with basemap.
        
All these routines were written 'on the fly', without a design plan
or test cases. Everything probably could use a rewrite!
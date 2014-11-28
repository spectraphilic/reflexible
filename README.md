# reflexible #

Github: https://github.com/spectraphilic/reflexible

[![Travis CI](https://travis-ci.org/spectraphilic/reflexible.png?branch=master)](https://travis-ci.org/spectraphilic/reflexible)

## About ##

This module is being created to work with FLEXPART model
ouput. FLEXPART is a Lagrangian particle dispersion model written by
Andreas Stohl. It is in active development at NILU and by Sabine
Eckhardt, John Burkhart and numerous other research groups worldwide.

This module is based partially on Matlab(c) routines developed
by Sabine Eckhardt and Harald Sodemann at NILU. Subsequent development
and the Python routines have been conducted by John Burkhart.

This software was written for the purposes of conducting research and displaying
model results for publications. More information can be found on my
[homepage](http://niflheim.nilu.no/~burkhart) under the 'Software' and
'Publications' sections.

## Contact, Issues, and Discussion ##

A maillist has been set up for users at sourceforge. Please direct
inquiries via the
[reflexible_maillist](https://lists.sourceforge.net/lists/listinfo/reflexible-users)

As for Issues, please use the
q[tool](https://github.com/spectraphilic/reflexible/issues) at
github.

Direct inquiries may be sent to::

	John F Burkhart, University of Oslo
	[john.burkhart](mailto://john.burkhart@geo.uio.no)

## Web Resources ##

### BitBucket ###

Clone the github repository into a folder in your PYTHONPATH::

  $ git clone https://github.com/spectraphilic/reflexible.git

### Sphinx Documentation ###

The documentation is hosted at
[readthedocs](http://reflexible.readthedocs.org/en/latest/index.html).

## Working with reflexible

There are a few 'gotchas' when using the module. First, you will
likely have to recompile (f2py) the FortFlex.f file and create a
FortFlex.so module for whatever computer you're using. For this you
need a Fortran compiler (gfortran) installed, so sometimes that can be
an issue.

An alternative 'BinaryFile' class has been created so one can work
with pure Python. Alone, it is significantly slower than the FortFlex
module, however, if you use the dumpgrid module, significant speedups
can be achieved. For this run::

  $ python setup.py build_ext --inplace

This will compile the pflexcy.pyx file into a pflexcy.so module that
can be imported and used by the pf.readgrid function. A series of
tests are run to try and determine the best module to use -- somewhat
transparently to the user. See the reflexible.readgridV8 function for
more information.

Primary functionality comes from the readheaderV8 function and the
readgridV8 function. I have created a "Header" class that can be used
for some typical analysis. See the examples directory, and don't
forget to read the source code and the [getting
started](http://reflexible.readthedocs.org/en/latest/getting_started.html)
documentation.

## Testing reflexible

reflexible comes with a suite of test units that you can run in a
series of ways after you compiled the extensions and before you
install it::

  $ PYTHONPATH=. python -c "import reflexible; reflexible.test()"

or::

  $ PYTHONPATH=. python reflexible/tests/all.py

or using the excellent py.test (recommended for developers)::

  $ PYTHONPATH=. py.test

## Installation

If all is working correctly, then you can install it (you might need
to be superuser here)::

  $ python setup.py install

And *hopefully* everything works!.

## Reporting problems

When you run into problems it is always nice that when you are filing
a ticket you would add the information about the versions you are
using.  You can do that via the `pf.print_versions()`.  Here it is an
output example::

  $ python -c "import reflexible; reflexible.print_versions()"
  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  reflexible version: 0.10.0
  NumPy version:     1.9.0
  Python version:    2.7.8 |Anaconda 2.1.0 (64-bit)| (default, Aug 21 2014, 18:22:21)
  [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]
  Platform:          linux2-x86_64
  Byte-ordering:     little
  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

## Key Tools

The primary workflow and usefulness of reflexible comes from the
read_header and read_grid routines, which are designed to help
ultimately with plotting and data analysis of the FLEXPART output.

To get started, at least become familiar with:

* Header class
* read_header
* read_grid
* fill_backward
* plot_sensitivity

NOTE: the other plotting functions are mostly wrappers which
ultimately make a call to plot_sensitivity.

## Disclaimer

About the mapping.py module, it is primarily designed as some
convenience routines for working with basemap.

All these routines were written 'on the fly', without a design plan or
test cases. Everything probably could use a rewrite!

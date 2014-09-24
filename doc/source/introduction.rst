.. _introduction:


************
Introduction
************

pflexible is an open source Python package to work with Lagrangian Particle
Disperion Model output. Currently it is built for `FLEXPART
<http://transport.nilu.no/flexpart>`_ but future versions will include greater
generality.

Contributions and collaboration are welcome. The code is hosted at bitbucket and
the documentation is hosted at readthedocs. pflexible is licensed under
Creative Commons.

Current development activities are focused on improved generality and handling of
FLEXPART output in all possible run configurations, with and without deposition,
forward, backward, or otherwise.

.. _getting-pflexible:

Getting pflexible
=================
Please contact: `John F. Burkhart <jfburkhart@gmail.com>`_

First, make sure you also have the dependencies installed:
  - numpy
  - matplotlib
  - basemap (matplotlib toolkit)
  - f2py (to build FortFlex)
  - netCDF4 (*not* python-netcdf)
  - PIL

Note the easiest way I've found to deal with the dependencies is to use one of
'complete distributions' such as `Enthought <http://www.enthought.com>`_ or the
`python(xy) <http://www.pythonxy.com>`_ or ideally `Anaconda <http://www.continuum.io>`_.
For Ubuntu you can pretty easily just install the required packages. For the netcdf, we 
found it easiest to use the science-meteorology-dev meta package and use pip to install netcdf4-python.

Once you've installed all the dependencies, you can get the code from either
sources below.

BitBucket
---------

The code is available to the public at `bitbucket
<https://bitbucket.org/jfburkhart/pflexible>`_. 

PyPi
----

The `pflexible <http://pypi.python.org/pypi/pflexible>`_ code is also posted to
pypi, but this is more likely to fall out of date.


----

.. _email-list:

Mailing List
============

There is a mailing list for the project set up at with sourceforge. You can subscribe to
the `pflexible
<https://lists.sourceforge.net/lists/listinfo/pflexible-users>`_ list for user
discussions.



.. _setting-pythonpath:

Installation and setting the PYTHONPATH
=======================================

If all is working correctly, and you have all the required dependencies, the it
should simply be a matter of running setup.py::

    python setup.py install


Depending on where you checked out the pflexible module to, you need to make
sure it is accessbile in your PYTHONPATH environment variable. The dependencies
also need to be available in the paths defined here. 

.. note::
   This can be accomplished after you've checked out the software::
      %export PYTHONPATH=/path/to/pflexible

.. sidebar:: Python at NILU

   Setting paths for most dependencies can be accomplished with::
   
    %export PYTHONPATH=$PYTHONPATH:/xnilu_wrk/jfb/hg

   And *hopefully* everything works!.

.. _building-FortFlex:

----

Building FortFlex
=================
 
FortFlex is a Fortran module that allows highly efficient reading of the raw FLEXPART output.
Building FortFlex is simple, and required only running a script within the ``f2py_build`` directory of 
the pflexible distribution::

    cd f2py_build
    chmod +x build_FortFlex.sh
    ./build_FortFlex.sh
    
Assuming you have all dependencies for f2py (gfortran, etc.) this will build and copy a new module into
pflexible package called FortFlex.so that will automatically be used by the ``pf.read_grid`` routines.

If you have built the ``FortFlex.so`` module independently, be sure to copy it into the same directory as
``pflexible.py`` or somewhere on your *PYTHONPATH*. **NOTE** I am trying to replace this dependency 
(or at least make more builds available, if you have suggestions, please contact me).


.. toctree::

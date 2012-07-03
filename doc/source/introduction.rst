.. _introduction:


***************
Introduction
***************

pflexible is an open source Python package to work with Lagrangian Particle
Disperion Model output. Currently it is built for `FLEXPART
<http://transport.nilu.no/flexpart>`_ but future versions will include greater
generality.

Contributions and collaboration are welcome. The code is hosted at bitbucket and
the documentation is hosted at readthedocs. pflexible was written by `John
F. Burkhart <http://niflheim.nilu.no/jfburkhart>`_ and is licensed under
Creative Commons.

.. _getting-pflexible:

Getting pflexible
=================
Please contact: `John F. Burkhart <jfburkhart@gmail.com>`_

The code is available to the public, but it is presently maintained at bitbucket
with restricted sharing. Upon request, we can bundle the latest release.

Make sure you also have the dependencies installed:
  - numpy
  - matplotlib
  - basemap (matplotlib toolkit)
  - f2py (to build FortFlex)
  - Scientific (not for long)
  - netCDF4
  - PIL

----
.. _setting-pythonpath:

Installation and setting the PYTHONPATH
-------------------------------------------------

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

Unless you're working on a 64bit Ubuntu machine, you will certainly have to
rebuild the FortFlex module. This is fairly easy. Navigate to the ``f2py_build``
directory. See the README file for instructions, but most likely, you just need
to run the following command::

   %f2py -c --fcompiler=gfortran FortFlex.pyf FortFlex.f

You will need to replace the `fcompiler` with whatever compiler you have. Once
you have built the ``FortFlex.so`` file, copy into the same directory as
``pflexible.py`` or somewhere on your *PYTHONPATH*. **NOTE** I am
trying to replace this dependency (or at least make more builds available, if
you have suggestions, please contact me).


.. toctree::

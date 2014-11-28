.. _introduction:


************
Introduction
************

reflexible is an open source Python package to work with Lagrangian
Particle Disperion Model output. Currently it is built for `FLEXPART
<http://transport.nilu.no/flexpart>`_ but future versions will include
greater generality.

Contributions and collaboration are welcome. The code is hosted at
github and the documentation is hosted at readthedocs. reflexible is
licensed under Creative Commons.

Current development activities are focused on improved generality and
handling of FLEXPART output in all possible run configurations, with
and without deposition, forward, backward, or otherwise.

.. _getting-reflexible:

Getting reflexible
=================
Please contact: `John F. Burkhart <jfburkhart@gmail.com>`_

First, make sure you also have the dependencies installed:
  - numpy
  - matplotlib
  - basemap (matplotlib toolkit)
  - f2py (to build FortFlex)
  - netcdf4-python (*not* python-netcdf)
  - PIL

Note the easiest way I've found to deal with the dependencies is to
use one of 'complete distributions' such as `Enthought
<http://www.enthought.com>`_ or the `python(xy)
<http://www.pythonxy.com>`_ or ideally `Anaconda
<http://www.continuum.io>`_.  For Ubuntu you can pretty easily just
install the required packages. For the netcdf, we found it easiest to
use the science-meteorology-dev meta package and use pip to install
netcdf4-python.

Once you've installed all the dependencies, you can get the code from either
sources below.

GitHub
------

The code is available to the public at `github
<https://github.com/spectraphilic/reflexible>`_.

PyPi
----

The `reflexible <http://pypi.python.org/pypi/reflexible>`_ code is
also posted to pypi, but this is more likely to fall out of date.


----

.. _email-list:

Mailing List
============

There is a mailing list for the project set up at with
sourceforge. You can subscribe to the `reflexible
<https://lists.sourceforge.net/lists/listinfo/reflexible-users>`_ list
for user discussions.

.. _installation:

Installation
============

If all is working correctly, and you have all the required
dependencies, then it should simply be a matter of running setup.py::

    python setup.py build

and then install it (you might need to be superuser here)::

    python setup.py install

And *hopefully* everything works!.


.. toctree::

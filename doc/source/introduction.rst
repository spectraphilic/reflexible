************
Introduction
************

reflexible is an open source Python package to work with Lagrangian
Particle Dispersion Model output. Currently it is built for `FLEXPART
<http://transport.nilu.no/flexpart>`_ but future versions will include
greater generality.

Contributions and collaboration are welcome. The code is hosted at
github and the documentation is hosted at readthedocs. reflexible is
licensed under Creative Commons.

Current development activities are focused on improved generality and
handling of FLEXPART output in all possible run configurations, with
and without deposition, forward, backward, or otherwise.

Getting reflexible
=================

The code is available to the public at `github
<https://github.com/spectraphilic/reflexible>`_.  You can easily clone
the git repository::

    $ git clone https://github.com/spectraphilic/reflexible.git

Build and test
==============

It should simply be a matter of changing the repo directory and running
setup.py::

    $ python setup.py build_ext --inplace

and then run the tests with::

    $ pytest

If the test suite pass, you can proceed with installation.

Installation
============

It should simply be a matter of running::

    $ python setup.py install

Quick install
=============

You may also want to install the package in one single shot (no testing
though!) with::

    $ pip install git+https://github.com/spectraphilic/reflexible.git

And *hopefully* everything works!.


.. toctree::

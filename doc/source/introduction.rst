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

A brief description of the package
=================================

The reflexible package is developed to work with output from the Lagrangian
Particle Dispersion Model, `FLEXPART <https://www.flexpart.eu/>`_ .

The module relies extensively on the users knowledge of FLEXPART data in
general, and thus one is strongly encouraged to read the
`users guide <http://zardoz.nilu.no/~andreas/flexpart/flexpart8.pdf>`_
which explains some basics regarding the model.

*Note* If you are interested in contributing functionality for other FLEXPART
versions, please contact me at `jfburkhart@gmail.com`

Purpose
-------

The purpose of the module is to make the creation of some standard plotting
products as easy as possible. However, due to the complex nature of FLEXPART
output, this isn't so easy! Regardless, I hope you find some of the
functionality helpful. The most critical functions are readheader and readgrid
which will at least get the data into Python so you can play with it as you are
most comfortable.

.. warning::
    You are entering the domain of a scientist trying to write code. Constructive
    input is sought, but don't complain if something breaks!


Getting reflexible
=================

The code is available to the public at `github
<https://github.com/spectraphilic/reflexible>`_.  You can easily clone
the git repository::

    $ git clone https://github.com/spectraphilic/reflexible.git

Install the requirements:

    $ conda install --file requirements.txt -c conda-forge

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

**NOTE** It is only planned to support Python 3.

Quick install
=============

You may also want to install the package in one single shot (no testing
though!) with::

    $ pip install git+https://github.com/spectraphilic/reflexible.git

And *hopefully* everything works!.


.. toctree::

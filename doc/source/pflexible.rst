.. _pflexible-overview:


********************
The pflexible module
********************

.. _module-description:

A brief description of the module
=================================

The pflexible module is developed to work with output from the Lagrangian
Particle Dispersion Model, `FLEXPART <http://transport.nilu.no/flexpart>`_ .

The module relies extensively on the users knowledge of FLEXPART data in
general, and thus one is strongly encouraged to read the `users guide <http://zardoz.nilu.no/~andreas/flexpart/flexpart8.pdf>`_ which explains some basics regarding the model.

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


*****************
The pflexible API
*****************

  :Release: |version|
  :Date: |today|
  :Author: John F. Burkhart

.. automodule:: pflexible
   :members:
   :undoc-members:

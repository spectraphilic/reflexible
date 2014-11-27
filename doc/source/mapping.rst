
.. _mapping-overview:


******************
The mapping module
******************

.. _module-description:

A brief description of the module
=================================

The mapping module is a helper function to the :mod:`reflexible` module.
Primarily it is designed to perform a few tasks relating to using the 
matplotlib `Basemap
<http://matplotlib.sourceforge.net/basemap/doc/html/api/basemap_api.html#module-mpl_toolkits.basemap>`_ 
module. I haven't confirmed whether how I pass the figures around or not is a
good idea, and would welcome suggestions.

.. warning::
  This module is not fully prepared for public use. There are a lot of
  custom functions, not written in a generic sense. Use with caution.


Purpose
-------

The purpose of this module is to ease create some basic mapping routines using
the basemap module. These are called directly from the :mod:`reflexible` for
example in the :func:`plot_sensitivity` routine. The core idea is that a
"FIGURE" object is created using the :func:`get_FIGURE` function which has some
key attributes. In general, this is transparent to the user, just intialize
a FIG object as NONE, then pass it to the functions with the FIGURE argument
set to your 'FIG' object.::

  > FIG = None
  > FIG = mp.plot_function(data,FIGURE=FIG)
  > 

The 'FIG' object can then be passed around and reused saving time and
resources. In general, the FIGURE object has the following attributes:

===============         ==================================
attribute / key         description
===============         ==================================
fig                     A fig object, use
                        plt.figure(FIG.fig.number) to make
                        it active
m                       A basemap instance for the plot
ax                      The primary axis instance
indices                 See the :func:`get_FIGURE` which 
                        describes the indices. 
===============         ==================================

Regions
-------
Another commonly used paradigm is the passing of a 'map_region' keyword to the
functions. Regions are defined manually at present. You'll have to edit the
:file:`mapping.py` and specifically, the :func:`map_regions`. Following the
instructions for the `Basemap
<http://matplotlib.sourceforge.net/basemap/doc/html/api/basemap_api.html#module-mpl_toolkits.basemap>`_
toolkit you can define your own unique region. See other regions as examples.


**Warning**
-----------

This is a module in active development, and there are no guarantees for backward
compatability. Constructive input is sought, but don't complain if something breaks!

.. toctree::

.. _api-index:

***************
The mapping API
***************

  :Release: |version|
  :Date: |today|

.. automodule:: mapping
   :members:
   :undoc-members:

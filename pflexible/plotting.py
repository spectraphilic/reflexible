# ### FLEXPART Plotting Functions  ########

import pdb

import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
# Basemap
try:
    from mpl_toolkits.basemap import shiftgrid, addcyclic
except ImportError:
    from matplotlib.toolkits.basemap import shiftgrid, addcyclic

from .helpers import Structure


def _plot_dropm():
    try:
        del FIGURE.m
        plt.close('all')
        del m
        "m dropped"
        return 1
    except:
        print 'Cannot delete m, does not exist.'
        return 0


def plot_trajchar(H, R, numrelease=1, varindex=[4, 15, 14],
                  FIGURE=None, map_region=None, projection=None, coords=None):
    """ plot FLEXPART trajectory (or set of trajectories) characteristics
        R is a dictionary returned from pflexible.read_trajectories

        varindex is a list of the variables to plot from the trajectories
    """
    if FIGURE is None:
        FIGURE = mp.get_FIGURE(getm=False)

    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    try:
        del fig.axes[0]
    except:
        print 'ax not removed'

    T = R['Trajectories']
    labels = R['labels']
    t = T[np.where(T[:, 0] == numrelease), :][0]
    numplot = len(varindex)

    for i in range(numplot):
        sp = fig.add_subplot(numplot, 1, i + 1)
        sp.plot(t[:, 1], t[:, varindex[i]])
        ax = plt.gca()
        if i + 1 != numplot:
            plt.setp(ax, xticklabels=[])
        plt.ylabel(labels[varindex[i]])
        plt.grid('on')

    return FIGURE


def plot_releases(R, FIGURE=None, threedim=False,
                  map_region=None, projection=None, coords=None,
                  overlay=True,
                  draw_circles=True,
                  draw_labels=True,
                  cbar2=True,
                  MapPar=None):
    """ plot the llon,llat,elv1 of the releases.

    Usage::
        >F = plot_releases(R)

    See the :func:`read_releases` function for information on how to
    generate "R"
    """

    if threedim:
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except:
            raise ImportError('mplot3d not available from mpl_toolkits!')

        fig = plt.figure()
        ax = Axes3D(fig)

        ax.plot(R.lllon, R.lllat, R.elv1)

        return fig
    else:
        # Set up the FIGURE
        if FIGURE is None:
            FIGURE = mp.get_FIGURE(map_region=map_region,
                                   projection=projection,
                                   coords=coords, MapPar=MapPar)
        # Get fig info and make active
        fig = FIGURE.fig
        m = FIGURE.m
        ax = FIGURE.ax
        plt.figure(fig.number)
        plt.axes(ax)

        # prepare the data
        lon = R.lllon
        lat = R.lllat
        zlevel = R.elv1
        zsize = np.ones(len(lon)) * 30
        marker = 'o'

        # clear the previous track
        if 'circles' in FIGURE.keys():
            del FIGURE['circles']
        if overlay is False:
            del ax.collections[FIGURE.indices.collections:]
            del ax.texts[FIGURE.indices.texts:]

        # plot the track
        cx, cy = m(lon, lat)
        if draw_circles:
            cmap = plt.get_cmap('gist_gray')
            circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                                marker=marker, edgecolor=None,
                                zorder=10, alpha=0.85)
            # m.scatter has no color bar,
            # so create a ghost 'scatter' instance:
            pos = ax.get_position()
            l, b, w, h = getattr(pos, 'bounds', pos)
            jnkfig = plt.figure()
            jnkax = jnkfig.add_axes([l, b, w, h], frameon=False)
            jnkmap = jnkax.scatter(cx, cy, zsize, zlevel,
                                   cmap=cmap, marker=marker,
                                   edgecolor=None)

            try:
                FIGURE.circles = circles
            except:
                pass
            # make the figure active again
            plt.figure(fig.number)
            # draw the legend and title
            # check if a colorbar legend exists (this will cause
            # problems for subplot routines!)
            if cbar2 is True:
                cax = plt.axes([l + w + 0.12, b, 0.02, h - 0.035])
            else:
                cax = plt.axes([l + w + 0.03, b, 0.025, h - 0.035])
            cb = fig.colorbar(jnkmap, cax=cax)  # draw colorbar
            p_leg = mpl.font_manager.FontProperties(size='6')
            cax.set_title('altitude\n(m)', fontproperties=p_leg)

            # delete the ghost instance
            plt.close(jnkfig.number)
            del jnkax, jnkfig, jnkmap

        plt.axes(ax)
        plt.setp(ax, xticks=[], yticks=[])

        if draw_labels:
            p_cax = mpl.font_manager.FontProperties(size='8',
                                                    style='italic',
                                                    weight='bold',
                                                    )

        FIGURE.fig = fig
        FIGURE.m = m
        FIGURE.ax = ax
        return FIGURE


def plot_spectra(inspectra,
                 plt_title='', fig_title='',
                 y_label=None,
                 spectra_label='bins',
                 FIGURE=None, y_datarange=None,
                 cum=False, labels=[],
                 bars=False, debug=False):
    """ plot a spectra

    Usage::

        > FIG = plot_agespectra(spectra,*kwargs)


    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      runid                 an ID for the title
      yunits                units for the y-axis [None]
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      data_range            y-axis data range
      web_version           Set to True if using the make_webpages
                            function. see: :func:`read_agespectrum`
      continental           Set to True if continental_spectrum
      cum                   Set to true if data should be cumulative
      bars                  Plot using stacked bars rather than a
                            filled line plot
      =============         ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes,
        maybe there is a
        better way to handle this.


    """
    # make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    if FIGURE is None:
        FIGURE = Structure()
        fig = plt.figure(figsize=(8, 6))
        FIGURE.fig = fig
        ax = fig.add_subplot(111)  # ,pos=[0.1,0.2,.8,.7])
        FIGURE.ax = ax
    else:
        fig = FIGURE.fig
        ax = FIGURE.ax

    try:
        numageclasses = inspectra.numageclass
        releasetimes = inspectra.releasetimes
        inspectra = inspectra[:]
    except:
        numageclasses = inspectra.shape[1] - 1
        releasetimes = inspectra[:, 0]
        inspectra = inspectra[:, 1:]

    if cum == 'norm':
        # Normalize the data so it fills to 100%
        # spectra = np.zeros(inspectra.shape)
        spectra = (
            inspectra.transpose() / np.sum(inspectra, axis=1)).transpose()
        spectra = np.cumsum(spectra[:, :], axis=1)
        # sums = np.sum(inspectra,axis=1)
        # for i,elem in enumerate(inspectra):
        #    spectra[i,:] = elem/sums[i]
    elif cum is True and cum != 'norm':
        spectra = np.cumsum(inspectra[:, :], axis=1)
    else:
        spectra = inspectra

    # Set up plotting environment colors
    Nc = np.array([float(i) / numageclasses for i in range(numageclasses)])
    norm = mpl.colors.normalize(Nc.min(), Nc.max())
    # jet = plt.cm.get_cmap('jet')
    jet = _gen_flexpart_colormap()
    plt.hold('on')
    facecolors = []
    # for i in range(0,H.numageclasses-1):
    for i in range(numageclasses):
        facecolors.append(jet(norm(Nc[i])))
        if labels:
            lbl = labels[i]
        else:
            lbl = i
        # ax.plot(ageclass[:,i])
        if i == 0:
            # create the baseline
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1])
            else:
                ax.fill_between(releasetimes, np.zeros(len(spectra[:, i])),
                                spectra[:, i],
                                color=facecolors[-1], label='%s' % lbl)
        else:
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1],
                       bottom=spectra[:, i - 1])
            else:
                ax.fill_between(releasetimes, spectra[:, i - 1], spectra[:, i],
                                color=facecolors[-1], label='%s' % (lbl))
                # facecolors.append(jet(norm(Nc[i+1])))

    # ax.set_yscale('log')
    if y_datarange:
        print 'setting data range'
        ax.set_ylim(y_datarange)
    else:
        ax.set_ylim((0, spectra.max()))
    ax.grid(True)
    if y_label:
        ax.set_ylabel('%s' % y_label, rotation=90, ha='center', size='small')

    if plt_title:
        plt.title('%s' % (plt_title), size='small')
    if fig_title:
        fig.suptitle('Emissions by %s' % (spectra_type), size='small')

    fig.autofmt_xdate()
    # for xl in ax.get_xticklabels():
    #    plt.setp(xl,size='x-small')
    # ListedColormap
    pos = ax.get_position()
    l, b, w, h = getattr(pos, 'bounds', pos)
    # polygons = ax.collections
    # for i,p in enumerate(polygons): p.set_color(facecolors[i])
    # BoundaryNorm, and extended ends to show the "over" and "under"
    # value co
    ax2 = fig.add_axes([l, .08, w, 0.03])
    cmap = mpl.colors.ListedColormap(facecolors)
    # cmap.set_over('0.25')
    # cmap.set_under('0.75')

    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.
    bounds = range(numageclasses + 1)
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
                                    #                                norm=norm,
                                    # to
                                    # use 'extend', you must
                                    #                                #
                                    # specify two extra boundaries:
                                    boundaries=bounds,
                                    #
                                    # extend='both',
                                    # values=range(H.numageclasses+1),
                                    ticks=bounds,  # optional
                                    spacing='proportional',
                                    orientation='horizontal')
    if spectra_label:
        cb2.set_label(spectra_label, size='x-small')
    if labels:
        ax2.set_xticklabels(labels, va='center', ha='left', rotation=0,
                            size='xx-small')

    if debug:
        plt.show()

    # make main axes current
    fig.add_subplot(ax)

    FIGURE.ax = ax
    FIGURE.fig = fig

    return FIGURE


def plot_agespectra(H, agespectra,
                    plt_title='', y_label=None,
                    cbar_labels=None,
                    FIGURE=None, y_datarange=None,
                    web_version=False, continental=False, cum=False,
                    bars=False, debug=False):
    """ plot an agespectra

    Usage::

        > FIG = plot_agespectra(H,agespectra,*kwargs)

    This creates an filled in between agespectra plot using information from
    the header "H".

    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      runid                 an ID for the title
      yunits                units for the y-axis [None]
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      data_range            y-axis data range
      web_version           Set to True if using the make_webpages
                            function. see: :func:`read_agespectrum`
      continental           Set to True if continental_spectrum
      cum                   Set to true if data should be cumulative
      bars                  Plot using stacked bars rather than a
                            filled line plot
      =============         ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes,
        maybe there is a
        better way to handle this.

    .. note::
        Required attributes of 'H' if you want to make a dummy H. Or just set
        H to "None" and it will be taken from the agespectra input.
        H.numageclasses
        H.releasetimes


    """
    # make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    if FIGURE is None:
        FIGURE = Structure()
        fig = plt.figure()
        FIGURE.fig = fig
        ax = fig.add_subplot(111)
        FIGURE.ax = ax
    else:
        fig = FIGURE.fig
        ax = FIGURE.ax

    if web_version:
        inspectra = agespectra.agespectrum[:, 4:]
    else:
        inspectra = agespectra[:]

    if continental:
        from nilu.pflexpart import emissions as em

        spectra_label = "Continents"
        spectra_type = "Continent"
        conts = em.Continents()
        cbar_labels = [i[1] for i in conts.flexpart_continents]
    else:
        spectra_label = "Ageclasses (Days)"
        spectra_type = "Ageclass"

    if not H:
        H = Structure()
        try:
            numageclasses = inspectra.numageclass
            releasetimes = inspectra.agespectrum[:, 0]
        except:
            numageclasses = inspectra.shape[1] - 1
            releasetimes = inspectra[:, 0]
            inspectra = inspectra[:, 1:]
    else:
        numageclasses = H.numageclasses
        releasetimes = H.releasetimes

    # if cum == 'norm':
    #    ## Normalize the data so it fills to 100%
    #    #spectra = np.zeros(inspectra.shape)
    #    spectra = (inspectra.transpose()/np.sum(inspectra, axis=1)).transpose()
    #    spectra = np.cumsum(spectra[:,:],axis=1)
    #    #sums = np.sum(inspectra,axis=1)
    #    #for i,elem in enumerate(inspectra):
    #    #    spectra[i,:] = elem/sums[i]
    # elif cum is True and cum != 'norm':
    #    spectra = np.cumsum(inspectra[:,:],axis=1)
    # else:
    #    spectra = inspectra
    if cum is not None:
        spectra = _cum_spec(inspectra, cum=cum)

    # Set up plotting environment colors
    Nc = np.array([float(i) / numageclasses for i in range(numageclasses)])
    norm = mpl.colors.normalize(Nc.min(), Nc.max())
    # jet = plt.cm.get_cmap('jet')
    jet = _gen_flexpart_colormap()
    plt.hold('on')
    facecolors = []
    # for i in range(0,H.numageclasses-1):
    for i in range(numageclasses):
        facecolors.append(jet(norm(Nc[i])))
        # ax.plot(ageclass[:,i])
        if i == 0:
            # create the baseline
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1])
            else:
                ax.fill_between(releasetimes, np.zeros(len(spectra[:, i])),
                                spectra[:, i],
                                color=facecolors[-1], label='%s' % i)
        else:
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1],
                       bottom=spectra[:, i - 1])
            else:
                ax.fill_between(releasetimes, spectra[:, i - 1], spectra[:, i],
                                color=facecolors[-1], label='%s' % (i))
                # facecolors.append(jet(norm(Nc[i+1])))

    # ax.set_yscale('log')
    if y_datarange:
        # print 'setting data range'
        ax.set_ylim(y_range)
    else:
        ax.set_ylim((0, spectra.max()))
    ax.grid(True)
    if y_label:
        ax.set_ylabel('%s' % y_label, rotation=90, ha='center', size='small')

    plt.title('%s' % (plt_title), size='small')
    fig.suptitle('Emissions by %s' % (spectra_type), size='small')
    fig.autofmt_xdate()
    # for xl in ax.get_xticklabels():
    #    plt.setp(xl,size='x-small')
    # ListedColormap
    pos = ax.get_position()
    l, b, w, h = getattr(pos, 'bounds', pos)
    # polygons = ax.collections
    # for i,p in enumerate(polygons): p.set_color(facecolors[i])
    # BoundaryNorm, and extended ends to show the "over" and "under"
    # value co
    ax2 = fig.add_axes([l, .08, w, 0.03])
    cmap = mpl.colors.ListedColormap(facecolors)
    # cmap.set_over('0.25')
    # cmap.set_under('0.75')

    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.
    bounds = range(numageclasses + 1)
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
                                    # norm=norm,
                                    # to use 'extend', you must
                                    # specify two extra boundaries:
                                    boundaries=bounds,
                                    # extend='both',
                                    # values=range(H.numageclasses+1),
                                    ticks=bounds,  # optional
                                    spacing='proportional',
                                    orientation='horizontal')
    cb2.set_label(spectra_label, size='x-small')
    if cbar_labels:
        cb2.set_ticklabels(cbar_labels)

    if debug:
        plt.show()

    return FIGURE


def plot_clusters(H, T, rel_i=0,
                  ncluster=0, sizescale=10000,
                  map_region=None, projection=None, coords=None,
                  FIGURE=None, overlay=False, opacity=0.7,
                  draw_circles=True, draw_labels=True,
                  MapPar=None):
    """ plots the clusters """
    trjs = T['Trajectories']
    labels = T['labels']

    # Set legend properties:
    p_legend = mpl.font_manager.FontProperties(size='8')

    # extract only releases of interest
    rel += 1  # account for zero indexing
    t = trjs[np.where(trjs[:, 0] == rel), :][0]
    if FIGURE is None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)

    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    # Remove prior retroplume
    if not overlay:
        if len(ax.artists) > 1:
            try:
                art_ind = len(ax.artists)
                del ax.artists[:]
                # del previous texts, but not lat/lon labels
                del ax.texts[-art_ind:]
            except:
                'could not delete artists'
                pass

    # get cluster from data
    # set up cluster index, start column of T array
    cindx = [16, 20, 24, 28, 32, 36]
    if ncluster == 'all':
        # plot all clusters
        ncluster = range(5)
    elif isinstance(ncluster, int):
        ncluster = [ncluster]

    for nc in ncluster:
        # clstrindx = 16 + (5*(nc))
        clstrindx = cindx[nc]
        indx = [1] + range(clstrindx, clstrindx + 4)
        data = t[:, indx]
        # use ellipses
        ells = _genEllipse(data, m, sizescale=sizescale)
        texts = []
        for item in ells:
            e = item[0]
            xy = item[1]
            texts.append([xy, e.get_label()])
            e.set_clip_box(ax.bbox)
            e.set_alpha(opacity)
            if draw_circles:
                ax.add_artist(e)
        for txt in texts:
            x = txt[0][0]
            y = txt[0][1]
            # print x,y, t[1]
            if draw_labels:
                ax.text(x, y, txt[1])
                # plt.colorbar(ax)
    # ax.legend(ax.artists,(o.get_label() for o in ax.artists), prop=p_legend )
    FIGURE.ax = ax
    FIGURE.fig = fig
    FIGURE.m = m

    return FIGURE


def plot_trajectory_ellipses(H, T, rel_i=0,
                             ncluster=0, sizescale=10000,
                             map_region=None, projection=None, coords=None,
                             FIGURE=None, overlay=False, opacity=0.7,
                             draw_circles=True, draw_labels=True,
                             MapPar=None):
    """
    Plots trajectories from FLEXPART output.

    """

    trjs = T['Trajectories']
    labels = T['labels']

    # Set legend properties:
    p_legend = mpl.font_manager.FontProperties(size='8')

    # extract only releases of interest according to rel_i
    rel = rel_i + 1  # account for zero indexing
    t = trjs[np.where(trjs[:, 0] == rel), :][0]
    if FIGURE is None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)

    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    # Remove prior retroplume
    if overlay:
        if len(ax.artists) > 1:
            try:
                art_ind = len(ax.artists)
                del ax.artists[:]
                # del previous texts, but not lat/lon labels
                del ax.texts[-art_ind:]
            except:
                'could not delete artists'
                pass

    # get trajectory from data
    indx = [1, 2, 3, 4]  # return, j,x,y,z
    data = t[:, indx]
    # use ellipses
    ells = _genEllipse(data, m, sizescale=sizescale)
    texts = []
    for item in ells:
        e = item[0]
        xy = item[1]
        texts.append([xy, e.get_label()])
        e.set_clip_box(ax.bbox)
        e.set_alpha(opacity)
        if draw_circles:
            ax.add_artist(e)
    for txt in texts:
        x = txt[0][0]
        y = txt[0][1]
        # print x,y, t[1]
        if draw_labels:
            ax.text(x, y, txt[1])
    plt.draw()
    FIGURE.ax = ax
    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.circles = ells

    return FIGURE


def plot_markers(H, lon, lat, zsize=None, zlevel=None,
                 FIGURE=None,
                 map_region=None, projection=None, coords=None,
                 overlay=True,
                 draw_circles=True,
                 draw_labels=False,
                 cbar2=True,
                 cbar2_title=None,
                 MapPar=None,
                 color='blue', edgecolor='none',
                 zorder=10, alpha=0.85):
    """Plot a group of x,y pairs, with optional size and color information.

    Usage::

        > FIG = plot_trajectory(H,RT,rel_id,*kwargs)

    .. note::
        You can substitude "None" for the :class:`Header` instance if you
        haven't
        created one


    Returns
      A :mod:`mapping` "FIGURE" object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         =============================================
      keyword               Description [default]
      =============         =============================================
      rel_i                 **required** release index
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      projection            A projection pre-defined in :mod:`mapping`
      coords                A set of lon,lat coords for autosetting
                            map map_region (not really working).
      overlay               [True] will overlay the trajectory
                            on top of another map instance.
      draw_circles          [True] will mark the trajectory with
                            circles. If [False] nothing is drawn.
      draw_labels           [False] unimplemented
      cbar2                 [True] draws the scale bar as a second axis.
      cbar2_title            Optional argument to overide the cbar title.
      MapPar                A Structure of mapping parameters to pass
                            to the basemap instance if desired.
      =============         =============================================

    .. todo::
        Who knows?! Could probably handle the overlaying better.

    .. note::
        Just set H to "None" if you haven't already created a "Header"
        instance.


    """

    if H:
        pass

    # Set up the FIGURE
    if FIGURE is None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)
    # Get fig info and make active
    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.fig.axes[0]
    plt.figure(fig.number)
    plt.axes(ax)

    # prepare the data
    cx, cy = m(lon, lat)
    if zlevel is None:
        zlevel = np.ones(len(lon)) * 10.
    if zsize is None:
        zsize = np.ones(len(lon)) * 1
    marker = 'o'

    # clear the previous track
    if 'circles' in FIGURE.keys():
        del FIGURE['circles']
    if overlay is False:
        del ax.collections[FIGURE.indices.collections:]
        del ax.texts[FIGURE.indices.texts:]

    # plot the track
    if draw_circles:
        cmap = plt.get_cmap('jet')
        circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                            marker=marker, edgecolor=edgecolor,
                            zorder=zorder, alpha=alpha)

        try:
            FIGURE.circles = circles
        except:
            pass
        # make the figure active again
        plt.figure(fig.number)
        # draw the legend and title
        # CREATE COLORBAR
        # does a colorbar already exist?
        # Get the current axes, and properties for use later
        ax0 = fig.axes[0]
        pos = ax0.get_position()
        l, b, w, h = pos.bounds
        if cbar2:
            try:
                cb2 = FIGURE.cb2
                cax2 = FIGURE.cax2
                cb2.update_normal(circles)
            except:
                cax2 = plt.axes([l + w + 0.08, b, 0.02, h])
                # using im2, not im (hack to prevent colors from being
                # too compressed at the low end on the colorbar - results
                # from highly nonuniform colormap)
                cb2 = fig.colorbar(circles,
                                   cax=cax2)  # , format='%3.2g') # draw
                # colorbar
                FIGURE.cax2 = cax2
                FIGURE.cb2 = cb2
                # check if a colorbar legend exists (this will cause
                # problems for subplot routines!)
            p_leg = mpl.font_manager.FontProperties(size='6')
            if cbar2_title:
                cax2.set_title(cbar2_title, fontproperties=p_leg)
            else:
                cax2.set_title('altitude\n(m)', fontproperties=p_leg)
        else:
            pass
            # cax2 = plt.axes([l+w+0.03, b, 0.025, h-0.2])
            # cb2 = fig.colorbar(jnkmap,cax=cax2) # draw colorbar

            # delete the ghost instance
            # plt.close(jnkfig.number)
            # del jnkax, jnkfig, jnkmap

    plt.axes(ax)
    plt.setp(ax, xticks=[], yticks=[])

    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.ax = ax
    try:
        FIGURE.cax2 = cax2
        FIGURE.cb2 = cb2
    except:
        pass
    return FIGURE


def plot_trajectory(H, T, rel_i, FIGURE=None,
                    map_region=None, projection=None, coords=None,
                    overlay=True,
                    draw_circles=True,
                    draw_labels=True, days_back=20,
                    cbar2=True,
                    cbar2_title=None,
                    MapPar=None):
    """Plot the center trajectory of the plume on the map

    Usage::

        > FIG = plot_trajectory(H,RT,rel_id,*kwargs)

    .. note::
        You can substitude "None" for the :class:`Header` instance if you
        haven't
        created one


    Returns
      A :mod:`mapping` "FIGURE" object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         =============================================
      keyword               Description [default]
      =============         =============================================
      rel_i                 **required** release index
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      projection            A projection pre-defined in :mod:`mapping`
      coords                A set of lon,lat coords for autosetting
                            map map_region (not really working).
      overlay               [True] will overlay the trajectory
                            on top of another map instance.
      draw_circles          [True] will mark the trajectory with
                            circles. If [False] nothing is drawn.
      draw_labels           [True] will draw the 'day' numbers on the
                            trajectories.
      days_back             For how many days back should the labels be
                            shown? [20]
      cbar2                 [True] draws the scale bar as a second axis.
      cbar2_title            Optional argument to overide the cbar title.
      MapPar                A Structure of mapping parameters to pass
                            to the basemap instance if desired.
      =============         =============================================

    .. todo::
        Who knows?! Could probably handle the overlaying better.

    .. note::
        Just set H to "None" if you haven't already created a "Header"
        instance.


    """

    if H:
        pass

    # Set up the FIGURE
    if FIGURE is None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)
    # Get fig info and make active
    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.fig.axes[0]
    plt.figure(fig.number)
    plt.axes(ax)

    # prepare the data
    trjs = T['Trajectories']
    rel = rel_i + 1  # # account for zero indexing

    # extract only releases of interest
    t = trjs[np.where(trjs[:, 0] == rel), :][0]

    # Get the data for the days_back we're interested in
    day_labels = _gen_daylabels(t[:days_back, 1])
    lon = t[:days_back, 2]
    lat = t[:days_back, 3]
    zlevel = t[:days_back, 4]
    zsize = np.ones(len(lon)) * 50
    marker = 'o'

    # clear the previous track
    if 'circles' in FIGURE.keys():
        del FIGURE['circles']
    if overlay is False:
        del ax.collections[FIGURE.indices.collections:]
        del ax.texts[FIGURE.indices.texts:]

    # plot the track
    cx, cy = m(lon, lat)
    if draw_circles:
        cmap = plt.get_cmap('gist_gray')
        circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                            marker=marker, edgecolor=None,
                            zorder=10, alpha=0.85)

        try:
            FIGURE.circles = circles
        except:
            pass
        # make the figure active again
        plt.figure(fig.number)
        # draw the legend and title
        # CREATE COLORBAR
        # does a colorbar already exist?
        # Get the current axes, and properties for use later
        ax0 = fig.axes[1]
        pos = ax0.get_position()
        l, b, w, h = pos.bounds
        if cbar2:
            try:
                cb2 = FIGURE.cb2
                cax2 = FIGURE.cax2
                cb2.update_normal(circles)
            except:
                cax2 = plt.axes([l + w + 0.08, b, 0.02, h])
                # using im2, not im (hack to prevent colors from being
                # too compressed at the low end on the colorbar - results
                # from highly nonuniform colormap)
                cb2 = fig.colorbar(circles,
                                   cax=cax2)  # , format='%3.2g') # draw
                # colorbar
                FIGURE.cax2 = cax2
                FIGURE.cb2 = cb2
        # check if a colorbar legend exists (this will cause
        # problems for subplot routines!)
        else:
            cax2 = plt.axes([l + w + 0.03, b, 0.025, h - 0.2])
            cb2 = fig.colorbar(jnkmap, cax=cax2)  # draw colorbar

        p_leg = mpl.font_manager.FontProperties(size='6')
        if cbar2_title:
            cax.set_title(cbar2_title, fontproperties=p_leg)
        else:
            cax2.set_title('altitude\n(m)', fontproperties=p_leg)
            # delete the ghost instance
            # plt.close(jnkfig.number)
            # del jnkax, jnkfig, jnkmap

    plt.axes(ax)
    plt.setp(ax, xticks=[], yticks=[])

    if draw_labels:
        p_cax = mpl.font_manager.FontProperties(size='10',
                                                style='italic',
                                                weight='bold',
                                                )

        for i, (p, x, y) in enumerate(zip(day_labels, cx, cy)):
            if x > m.llcrnrx and x < m.urcrnrx:
                if y > m.llcrnry and y < m.urcrnry:
                    ax.text(x, y + y * .02, '{0}'.format(p), va='bottom',
                            ha='left',
                            fontproperties=p_cax, zorder=11,
                            color='white',
                            bbox=dict(facecolor='green', alpha=0.5)
                            )

    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.ax = ax
    FIGURE.cax2 = cax2
    FIGURE.cb2 = cb2
    return FIGURE


def plot_at_level(H, D, level=1,
                  ID=' ',
                  map_region=5, projection='lcc',
                  overlay=False,
                  datainfo_str=None, log=True,
                  data_range=None, coords=None, FIGURE=None,
                  plot_title=None,
                  units=None,
                  **kwargs
                  ):
    """
    TODO:
    -make units a function of H['species']
    """
    if units is None:
        units = H.output_unit
    if D is None:
        D = H.D
    rel_i = D.rel_i
    species = D.species
    timestamp = D.timestamp
    data = D.slabs[level]

    if level == 1:
        level_desc = 'Footprint'
    else:
        level_desc = H.outheight[level - 1]

    dmax = data.max()
    dmin = data.min()
    # print dmax,dmin

    if data_range is None:
        data_range = [dmin, dmax]

    if H.direction == 'backward' and H.options['readp']:
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            # need to find a way to set m.a.s.l. or hPa here
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f,
            Z2: %.2f (%s)\n""" \
                           % (dmax, units, zp1, zp2, H.alt_unit)
        if plot_title is None:
            plot_title = """
        %s Sensitivity at %s %s: %s\n
        Release Start: %s, Release End: %s""" % \
                         (ID, level_desc, H['alt_unit'], species,
                          H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s """ % (dmax, units)
        if plot_title is None:
            plot_title = """ %s Sensitivity at %s %s: %s \n %s """ % (
                ID, level_desc, H['alt_unit'], species, timestamp)

    FIGURE = plot_sensitivity(H, data,
                              data_range=data_range,
                              rel_i=rel_i, log=log,
                              map_region=map_region, projection=projection,
                              units=units, datainfo_str=datainfo_str,
                              overlay=overlay,
                              coords=coords, FIGURE=FIGURE, **kwargs)

    FIGURE.ax.set_title(plot_title, fontsize=10)
    return FIGURE


plot_footprint = plot_at_level


def plot_totalcolumn(H, D=None,
                     ID=' ',
                     map_region=5, projection='lcc',
                     data_range=None, coords=None,
                     FIGURE=None, overlay=False,
                     datainfo_str=None, **kwargs):
    if D is None:
        D = H.D[(0, 0)]

    if 'units' in kwargs:
        units = kwargs.pop('units')
    else:
        units = 'ns m kg-1'

    rel_i = D.rel_i
    species = D.species
    timestamp = D.timestamp
    data = D.slabs[0]

    if data_range is None:
        dmax = data.max()
        dmin = data.min()
        data_range = [dmin, dmax]
    else:
        dmin, dmax, = data_range
        # print dmin,dmax

    if H.direction == 'backward':
        rel_i = D.rel_i
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f,
            Z2: %.2f (%s)\n""" % (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n
        Release Start: %s, Release End: %s""" % \
                     (ID, species, H['releasestart'][rel_i],
                      H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s""" % (dmax, units)
        plot_title = """
        %s Total Column Sensitivity: %s\n %s """ % (ID, species, timestamp)

    FIGURE = plot_sensitivity(H, data,
                              data_range=data_range,
                              rel_i=rel_i, map_region=map_region,
                              projection=projection, units=units,
                              datainfo_str=datainfo_str, coords=coords,
                              FIGURE=FIGURE, overlay=overlay, **kwargs)

    FIGURE.ax.set_title(plot_title, fontsize=10)

    return FIGURE


def plot_sourcecontribution(H, D,
                            rel_i=None, s=None,
                            ID=' ',
                            map_region='POLARCAT',
                            data_range=None, coords=None,
                            datainfo_str=None,
                            units=None,
                            FIGURE=None,
                            overlay=False,
                            cax_title=None,
                            **kwargs):
    """ plot_sourcecontribution: a wrapper around :func:`plot_sensitivity` for
    plotting a source contribution field.

    .. warning:
        This function is still quite finicky, I'm working on resolving it.

    Usage::

        > FIG = plot_sourcecontribution(H,D,*kwargs)

    Inputs
      H = a :class:`Header` instance from a FLEXPART run
      D = a grid instance from a FLEXDATA dictionary. See :func:`readgridV8`

      .. note::
          I am trying to create this so D can be a number of different types
          from a 2-d array to the grid instance type.


    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ================================================
      keyword               Description [default]
      =============         ================================================
      ID                    Can be passed to prefix the title.
      data_range            range of data for scale bars, if None will
                            be taken from min/max of data
      datainfo_str          A string for labeling the scale bar.
      cax_title             A string to pass to plot_sensitivity for colorbar
                            title (units will be passed as format argument)
      k                     Release index to plot from the data array (aka
      rel_i)
      s                     Species index
      map_region                A map_region specified in mapping.py
      projection            [deprecated] use pre-defined map_regions.
      overlay               Force removal of previous figure elements.
      FIGURE                A FIGURE instance from mapping module get_FIGURE
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.

    .. note::
        k and rel_i are used throughout pflexible, should be more consistent
        in future.


    """
    if D is None:
        D = H.D[(0, 0)]
        data = D.slabs[1]
    elif isinstance(D, list):
        data = D[0]
    elif isinstance(D, np.ndarray):
        data = D
    else:
        try:
            D = D[(s, rel_i)]
            data = D.slabs[0]  # take total column
        except:
            raise IOError("Unrecognized 'grid' type passed")

    if s is None:
        try:
            species = D.species
        except:
            species = 'unknown'
    else:
        species = H.species[s]

    try:
        timestamp = D.timestamp
    except:
        if rel_i is None:
            try:
                timestamp = D.timestamp
                rel_i = D.rel_i
            except:
                raise IOError("Either a release index (rel_i) must be passed, \
                           or the grid must have a timestamp attribute")
        else:
            timestamp = H.releasetimes[rel_i]

    if units is None:
        units = 'emission units'

    dmax = data.max()
    dmin = data.min()

    if data_range is None:
        data_range = [dmin, dmax]

    if H.direction == 'backward':
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f ,
            Z2: %.2f (%s)\n""" % \
                           (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n
        Release Start: %s, Release End: %s""" % \
                     (ID, species, H['releasestart'][rel_i],
                      H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f ,
            Z2: %.2f (%s)\n""" % \
                           (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n %s""" % (ID, species, timestamp)

    FIGURE = plot_sensitivity(H, data,
                              data_range=data_range,
                              rel_i=rel_i, map_region=map_region,
                              units=units,
                              datainfo_str=datainfo_str, coords=coords,
                              FIGURE=FIGURE,
                              cax_title=cax_title,
                              overlay=overlay, **kwargs)
    FIGURE.ax.set_title(plot_title, fontsize=10)
    return FIGURE


def plot_sensitivity(H, data,
                     data_range=None,
                     units='ns m^2 / kg',
                     datainfo_str=None,
                     plottitle=None,
                     rel_i=None,
                     map_region=None, projection=None,
                     dropm=None, coords=None,
                     overlay=False,
                     transform=True,
                     log=True,
                     FIGURE=None,
                     MapPar=None,
                     FigPar=None,
                     cax_title=None,
                     autofit=False, method='contourf', lsmask=False):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::

        > FIG = plot_sensitivity(H,data,*kwargs)

    This returns the FIGURE object, and plots the sensitivity from the data
    contained in the "D"
    array.

    Inputs
      H = a :class:`Header` instance for a FLEXPART run.
      data = a 2d data array containing the sensitivity values to plot,
      this can be extracted from a
      grid instance (see :func:`readgridV8` and :func:`get_slabs`)

    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ================================================
      keyword               Description [default]
      =============         ================================================
      data_range            range of data for scale bars, if None will
                            be taken from min/max of data
      cax_title             string to be passed as colorbar title (units will
                            be passed to the format argument)
      units                 units for the scale bar
      datainfo_str          A string for labeling the scale bar.
      plottitle             Title for the plot.
      rel_i                 Release index to plot from the data array
      map_region                A map_region specified in mapping.py
      projection            [deprecated] use pre-defined map_regions.
      dropm                 Force creation of a new basemap instance
      coords                Used with autofit option. An array of lat,lon
                            values for the mapping module to autofit a
                            basemap instance to.
      autofit               Try to generate map_region automagically (flakey)
      overlay               Force removal of previous figure elements.
      transform             For use with imshow method, if your data is not
                            in same coordinates as projection, try to transform
                            the data to the basemap projection.
      log                   Create a logarithmic color scale.
      FIGURE                A FIGURE instance from mapping module get_FIGURE
      MapPar                A Structure of paramters to be passed to the
                            basemap class when creating an instance.
      method                The method to use for plotting array data. May be
                            one of: [pcolormesh], imshow, or contourf
      lsmask                set to True to draw a grey landseamask [False]
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.
        lsmask = True seems to only work with certain projections (POLARCAT)

    .. note::
        This is the primary main function for creating plots of flexpart
        output.
        Most the other routines are simply wrappers to this function, passing
        arguments in with some predefined settings. For more information on the
        mechanics of this function, see the mapping.py module and the
        matplotlib
        basemap toolkit.


    """

    methods = ['imshow', 'pcolormesh', 'contourf', 'contour', 'None']
    assert method in methods, "method keyword must be one of: %s" % methods

    if FIGURE is None:
        if map_region is None:
            if MapPar is None:
                if autofit:
                    MapPar = _gen_MapPar_fromHeader(H)

        FIGURE = mp.get_FIGURE(map_region=map_region,
                               projection=projection, coords=coords,
                               MapPar=MapPar, FigPar=FigPar)
    else:
        if FIGURE.m is None:
            FIGURE = mp.get_FIGURE(fig=FIGURE.fig, ax=FIGURE.ax,
                                   map_region=map_region,
                                   projection=projection, coords=coords,
                                   MapPar=MapPar, FigPar=FigPar)

    if overlay is False:
        del FIGURE.ax.images[FIGURE.indices.images:]
        del FIGURE.ax.collections[FIGURE.indices.collections:]
        del FIGURE.ax.lines[FIGURE.indices.lines:]

    if dropm is not None:
        try:
            del m
            plt.close('all')
        except:
            print 'could not drop m'

    # make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax

    # make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    # set up transformations for the data array
    if method == 'imshow':
        if m.projection not in ['cyl', 'merc', 'mill']:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)),
                             H.dyout)[:-1]
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)),
                             H.dxout)[:-1]
            data = data[:-1, :-1]
        else:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)),
                             H.dyout)
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)),
                             H.dxout)

        # transform to nx x ny regularly spaced native projection grid
        if transform:
            dx = 2. * np.pi * m.rmajor / len(lons)
            nx = int((m.xmax - m.xmin) / dx) + 1
            ny = int((m.ymax - m.ymin) / dx) + 1
            if nx is 1:
                topodat = data
            else:
                topodat = m.transform_scalar(data, lons, lats, nx, ny)
        else:
            topodat = data

    if method != 'imshow':
        # Check to see if a cyclic wraparound is required
        lons = np.arange(H.outlon0, H.outlon0 + (H.dxout * H.numxgrid),
                         H.dxout)
        lats = np.arange(H.outlat0, H.outlat0 + (H.dyout * H.numygrid),
                         H.dyout)
        # if required add polar coordinates
        if 'npstere' in m.projection:
            if lats[-1] != 90.:
                npole = np.ones(len(lons)).T * data[0, :]
                npole = np.reshape(npole, (1, -1))
                data = np.vstack((npole, data))
                lats = np.hstack((lats, [lats[-1] + H.dyout]))

        if m.projection != 'merc':
            if lons[-1] - lons[0] < 360.:
                topodat, lons = addcyclic(data, lons)

        if m.projection == 'merc':
            topodat = data

    # get min/max range
    if data_range is not None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = data_range(data)

    if log:
        clevs = _log_clevs(dat_min, dat_max)

    else:
        clevs = [i for i in
                 np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    # draw land sea mask
    # m.fillcontinents(zorder=0)
    if lsmask:
        m.drawlsmask(ocean_color='grey', zorder=-10)

    # Plot Release Location if points were read
    if H.options['readp']:
        if rel_i:
            releaselocation = (H.xpoint[rel_i], H.ypoint[rel_i])
            xpt, ypt = m(releaselocation[0], releaselocation[1])
            # Remove prior location point
            try:
                del ax.lines[-1]
            except:
                pass
            location, = m.plot([xpt], [ypt], 'bx', linewidth=6, markersize=20,
                               zorder=1000)

    # Plot the footprint

    # Set up the IMAGE
    # cmapnames = ['jet', 'hsv', 'gist_ncar', 'gist_rainbow', 'cool', 'spectral']
    # colmap = plt.get_cmap('jet')
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    # Plotting METHODS (pcolormesh now default, imshow is smoother)
    # print topodat.max(), topodat.min(), topodat.shape
    if method == 'imshow':
        im = m.imshow(topodat, cmap=colmap, zorder=-1,
                      norm=mpl.colors.LogNorm(vmin=clevs[0],
                                              vmax=clevs[-1]))

    if method == 'pcolormesh':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.pcolormesh(nx, ny, topodat, cmap=colmap,
                          norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                  vmax=clevs[-1]))
    if method == 'contourf':
        # Trying some fancier scaling
        # cnts,bins = np.histogram(topodat,bins=100)
        # topodat = np.ma.masked_where(topodat< .05* np.average((0,
        # bins[1])),topodat)
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    if method == 'contour':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contour(nx, ny, topodat, cmap=colmap,
                       norm=mpl.colors.LogNorm(vmin=clevs[0],
                                               vmax=clevs[-1]))

    # Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    # CREATE COLORBAR
    # Note, with upgrades to matplotlib and basemap had to make some
    # changes here... no more 'ghost' axes
    # does a colorbar already exist?
    try:
        cb = FIGURE.cb
        cax = FIGURE.cax
        cb.update_normal(im)
    except:
        # make a copy of the image object, change
        # colormap to linear version of the precip colormap.
        # pdb.set_trace()
        # im2 = copy.copy(im)
        # im2.set_cmap(colmap)
        # create new axis for colorbar.
        h = 0.5 * h
        l = l + w + .03
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
        # using im2, not im (hack to prevent colors from being
        # too compressed at the low end on the colorbar - results
        # from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax)  # , format='%3.2g') # draw colorbar
        FIGURE.cax = cax
        FIGURE.cb = cb
    # cb.update_normal(im2)

    # set colorbar label and ticks
    # pdb.set_trace()
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10])  # #clevs, by 10 steps
    clabels.append(clevs[-1])  # # add the last label
    # cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])
    # fontproperties=p_cax)
    if H.direction == 'forward':
        cax.set_title('%s' % units,
                      fontproperties=p_cax)
    else:
        if cax_title:
            cax.set_title(cax_title.format(units), fontproperties=p_cax)
        else:
            cax.set_title('sensitivity\n({0})'.format(units),
                          fontproperties=p_cax)

    # make the original axes current again
    plt.axes(ax)

    # write text information block on plot
    # first try to remove prior text by accessing
    # the last text element added to the axes from
    # the prior iteration.
    # This is tricky when using together with plot_clusters...
    # need to figure out how to resolve the indexing
    # of what texts, collections, etc to delete, when iterating.
    try:
        del ax.texts[FIGURE.indices.texts:]
        del ax.artists[FIGURE.indices.artists:]
    except:
        pass
    if datainfo_str:
        plt.text(l, b + 1000,
                 datainfo_str,
                 fontsize=10,
                 bbox=dict(boxstyle="round",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           alpha=0.8
                           )
                 )

    FIGURE.ax = ax
    FIGURE.m = m
    FIGURE.fig = fig

    if plottitle is not None:
        # plt.title(plottitle,fontproperties=p_cax)
        # plt = plt
        FIGURE.ax.set_title(plottitle, fontsize=10)
    return FIGURE


def plot_curtain(H, data,
                 nx=None,
                 ny=None,
                 data_range=None,
                 units='ppbv',
                 datainfo_str=None,
                 asl=True,
                 plottitle=None,
                 log=True,
                 FIGURE=None,
                 cax_title=None,
                 method='contourf',
                 figPar=None):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::

        > FIG = plot_sensitivity(H,data,*kwargs)

    This returns the FIGURE object, and plots the sensitivity from the data
    contained in the "D"
    array.

    Inputs
      H = a :class:`Header` instance for a FLEXPART run.
      data = a 2d data array containing the sensitivity values to plot,
      this can be extracted from a
      grid instance (see :func:`readgridV8` and :func:`get_slabs`)

    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ================================================
      keyword               Description [default]
      =============         ================================================
      data_range            range of data for scale bars, if None will
                            be taken from min/max of data
      asl                   [True] plot the units in m.a.s.l
      cax_title             string to be passed as colorbar title (units will
                            be passed to the format argument)
      units                 units for the scale bar
      datainfo_str          A string for labeling the scale bar.
      plottitle             Title for the plot.
      rel_i                 Release index to plot from the data array
      map_region                A map_region specified in mapping.py
      projection            [deprecated] use pre-defined map_regions.
      dropm                 Force creation of a new basemap instance
      coords                Used with autofit option. An array of lat,lon
                            values for the mapping module to autofit a
                            basemap instance to.
      autofit               Try to generate map_region automagically (flakey)
      overlay               Force removal of previous figure elements.
      transform             For use with imshow method, if your data is not
                            in same coordinates as projection, try to transform
                            the data to the basemap projection.
      log                   Create a logarithmic color scale.
      FIGURE                A FIGURE instance from mapping module get_FIGURE
      MapPar                A Structure of paramters to be passed to the
                            basemap class when creating an instance.
      method                The method to use for plotting array data. May be
                            one of: [pcolormesh], imshow, or contourf
      lsmask                set to True to draw a grey landseamask [False]
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.
        lsmask = True seems to only work with certain projections (POLARCAT)

    .. note::
        This is the primary main function for creating plots of flexpart
        output.
        Most the other routines are simply wrappers to this function, passing
        arguments in with some predefined settings. For more information on the
        mechanics of this function, see the mapping.py module and the
        matplotlib
        basemap toolkit.


    """

    methods = ['imshow', 'pcolormesh', 'contourf', 'contour', 'None']
    assert method in methods, "method keyword must be one of: %s" % methods

    if FIGURE is None:
        FIGURE = Structure()

        fig = plt.figure(**figPar)
        ax = fig.add_subplot(111)

        FIGURE['fig'] = fig
        FIGURE['ax'] = ax

    fig = FIGURE.fig
    ax = FIGURE.ax

    # make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    # get min/max range
    if data_range is not None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = data_range(data)

    if log:
        clevs = _log_clevs(dat_min, dat_max)
    else:
        clevs = [i for i in
                 np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    # Set up the IMAGE
    # cmapnames = ['jet', 'hsv', 'gist_ncar', 'gist_rainbow', 'cool',
    # 'spectral']
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    topodat = data

    if method == 'imshow':
        im = plt.imshow(np.flipud(topodat), cmap=colmap, zorder=-1,
                        norm=mpl.colors.LogNorm(
                            vmin=clevs[0], vmax=clevs[-1]))

    if method == 'pcolormesh':
        im = plt.pcolormesh(nx, ny, topodat, cmap=colmap,
                            norm=mpl.colors.LogNorm(
                                vmin=clevs[0], vmax=clevs[-1]))
    if method == 'contourf':
        im = plt.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                          norm=mpl.colors.LogNorm(
                              vmin=clevs[0], vmax=clevs[-1]))

    if method == 'contour':
        im = plt.contour(nx, ny, topodat, cmap=colmap,
                         norm=mpl.colors.LogNorm(
                             vmin=clevs[0], vmax=clevs[-1]))

    # Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    # CREATE COLORBAR
    # Note, with upgrades to matplotlib and basemap had to make some
    # changes here... no more 'ghost' axes
    # does a colorbar already exist?
    try:
        cb = FIGURE.cb
        cax = FIGURE.cax
        cb.update_normal(im)
    except:
        # make a copy of the image object, change
        # colormap to linear version of the precip colormap.
        # create new axis for colorbar.
        h = 0.8 * h
        l = l + w + .02
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
        # using im2, not im (hack to prevent colors from being
        # too compressed at the low end on the colorbar - results
        # from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax)  # , format='%3.2g') # draw colorbar
        FIGURE.cax = cax
        FIGURE.cb = cb

    # set colorbar label and ticks
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10])  # #clevs, by 10 steps
    clabels.append(clevs[-1])  # # add the last label
    # cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])
    # fontproperties=p_cax)

    if cax_title:
        cax.set_title(cax_title.format(units), fontproperties=p_cax)
    else:
        cax.set_title('sensitivity\n({0})'.format(units),
                      fontproperties=p_cax)

    # make the original axes current again
    plt.axes(ax)
    plt.grid(True)
    FIGURE.ax = ax
    FIGURE.fig = fig

    if plottitle is not None:
        FIGURE.ax.set_title(plottitle, fontsize=10)

    return FIGURE


def plot_METDATA(METDATA, FIGURE, date=None, level=None):
    """ plots met data returned from :module:`pflexible.mapping.get_OpenDAP`

    """

    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax

    # make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    datelabels = METDATA['extracted_dates']
    slpin = METDATA['prmslmsl']['data'] * .01
    uin = METDATA['ugrdprs']['data']
    vin = METDATA['vgrdprs']['data']
    longitudes = METDATA['longitudes']
    latitudes = METDATA['latitudes']

    if date:
        cnt = datelabels.index(date)
    else:
        cnt = 0
    # FOR OPENDAP OVERLAY
    slp_ = slpin[cnt, :, :]
    u_ = uin[cnt, :, :]
    v_ = vin[cnt, :, :]

    # plot wind vectors on projection grid (looks better).
    # first, shift grid so it goes from -180 to 180 (instead of 0 to 360
    # in longitude).  Otherwise, interpolation is messed uplt.

    if longitudes[-1] - longitudes[0] < 360.:
        longarray = np.array(longitudes)
        slp, cyclon = addcyclic(slp_, longarray)
        u, cyclon = addcyclic(u_, longarray)
        v, cyclon = addcyclic(v_, longarray)
        ugrid, newlons = shiftgrid(180., u, cyclon, start=False)
        vgrid, newlons = shiftgrid(180., v, cyclon, start=False)
        slpgrid, newlons = shiftgrid(180., slp, cyclon, start=False)
    # transform vectors to projection grid.

    lons, lats = np.meshgrid(newlons, latitudes)
    # set desired contour levels.
    clevs = np.arange(960, 1061, 5)
    # compute native x,y coordinates of grid.
    x, y = m(lons, lats)
    # pos = FIG.ax.get_position()
    # l, b, w, h = pos.bounds

    CS = m.contour(x, y, slpgrid, clevs, linewidths=1, colors='k',
                   animated=True)
    # CS = FIG.m.contourf(x,y,slpgrid,clevs,cmap=plt.cm.RdBu_r,animated=True,alpha=0.7)
    # CS = FIG.m.contour(x,y,slp[0,:,:],clevs,linewidths=0.5,colors='k',animated=True)
    #            CS = FIG.m.contourf(x,y,slp[0,:,:],clevs,cmap=plt.cm.RdBu_r,animated=True,alpha=0.7)
    # plt.clabel(CS,inline=1,fontsize=10)
    # plot wind vectors over maplt.
    urot, vrot, xx, yy = m.transform_vector(
        ugrid, vgrid, newlons, latitudes, 51, 51,
        returnxy=True, masked=True)
    Q = m.quiver(xx, yy, urot, vrot, scale=500)
    # make quiver key.
    qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')

    FIGURE.ax = ax
    FIGURE.m = m
    FIGURE.fig = fig

    return FIGURE

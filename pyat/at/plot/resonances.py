"""AT plotting functions related to resonances."""

from __future__ import annotations

__all__ = ["farey_sequence", "plot_tune_diagram", "create_linepalette"]

import warnings
from fractions import Fraction

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np

# 2024jul31 oblanco at ALBA CELLS


def create_linepalette(
    linestyle: str | dict = None,
    linecolor: str = None,
    linewidth: int = None,
    addtolabel: str = None,
) -> dict[str, any]:
    """Create a line palette to plot resonance lines.

    Parameters:
        linestyle: If a dictionary is passed, it should contain
            {"normal": style1, "skew": style2}

            If 'dots' it uses dotted styles as linestyles. Equivalent to:
            {"normal": "dashdot", "skew": "dotted"}

            Default: {"normal": '-', "skew": '--'}
        linecolor: line color, e.g. "k". Default: custom values.
            See :py:func:`plot_tune_diagram`
        linewidth: line width. Default: custom values. See :py:func:`plot_tune_diagram`
        addtolabel: adds a string to the line label

    Returns:
        palette: dictionary containing the line properties for resonance plots.
    """
    # create default dictionary with line properties
    mypalettecolor1 = {
        1: "k",
        2: "b",
        3: "r",
        4: "g",
        5: "m",
        6: "c",
        7: "y",
        8: "darkcyan",
        9: "lightgreen",
        10: (0.1, 0.1, 0.1),
        11: (0.1, 0.1, 0.1),
        12: (0.2, 0.2, 0.2),
        13: (0.3, 0.3, 0.3),
        14: (0.4, 0.4, 0.4),
        15: (0.5, 0.5, 0.5),
    }
    mypalettewidth_main = {
        1: 4,
        2: 3,
        3: 2,
        4: 1,
        5: 1,
        6: 1,
        7: 1,
        8: 1,
        9: 1,
        10: 1,
        11: 1,
        12: 1,
        13: 1,
        14: 1,
        15: 1,
    }
    mypalettestyle_main = {"normal": "-", "skew": "--"}
    mypalettestyle_withdots = {"normal": "dashdot", "skew": "dotted"}
    # set the line style
    mypalettestyle = mypalettestyle_main
    if linestyle == "dots":
        mypalettestyle = mypalettestyle_withdots
    if isinstance(linestyle, dict):
        mypalettestyle = linestyle
    # set the line color
    mypalettecolor = mypalettecolor1
    if linecolor is not None:
        for lsize in range(1, 16):
            mypalettecolor[lsize] = linecolor
    # set the line width
    mylinewidth = mypalettewidth_main
    if isinstance(linewidth, int) and linewidth > 0:
        for lsize in range(1, 16):
            mylinewidth[lsize] = linewidth
    # add to label
    if addtolabel is None:
        addtolabel = ""
    # set the dictionary
    linepropdict = {}
    resontypes = ["normal", "skew"]
    for resontype in resontypes:
        linepropdict[resontype] = {}
    for resontype in resontypes:
        for idx in range(1, 16):
            propaux = {
                "color": mypalettecolor[idx],
                "linestyle": mypalettestyle[resontype],
                "linewidth": mylinewidth[idx],
                "label": str(idx) + resontype[0] + addtolabel,
            }
            linepropdict[resontype][idx] = propaux
    return linepropdict


def farey_sequence(nthorder: int, verbose: bool = False) -> tuple[list, list]:
    """
    Return the Farey sequence, and the resonance sequence of nth order.

    Parameters:
        nthorder: natural number bigger than 0.
        verbose: prints extra info. Default: False.

    Returns:
        fareyseqfloat: list of elements with the Farey sequence in Float format.
            See Eqs.(1,2,3) of [1].
        fareyseqfrac: list of elements with the Farey sequence in Fraction format.
            See Eqs.(1,2,3) of [1].

    Raises:
        ValueError: if given order is lower than 0, or window is zero.

    [1] R.Tomas. 'From Farey sequences to resonance diagrams.
            Phys.Rev.Acc.Beams 17, 014001 (2014)'
    """
    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    if nthorder < 1:
        raise ValueError("The order must be a positive integer larger than zero.")
    verboseprint(f"nthorder={nthorder}")

    afarey = 0
    bfarey = 1
    cfarey = 1
    dfarey = nthorder
    farey = [0, 1 / dfarey]
    fracfarey = [Fraction(0), Fraction(1, dfarey)]
    idx = 0
    while (farey[-1] < 1) and (idx < 100):
        idx += 1
        caux = np.floor((nthorder + bfarey) / dfarey) * cfarey - afarey
        daux = np.floor((nthorder + bfarey) / dfarey) * dfarey - bfarey
        afarey = cfarey
        bfarey = dfarey
        cfarey = int(caux)
        dfarey = int(daux)
        farey.append(cfarey / dfarey)
        fracfarey.append(Fraction(cfarey, dfarey))
    verboseprint(f"farey_float{nthorder}= {farey}")
    verboseprint(f"farey_frac{nthorder} = {fracfarey}")
    return farey, fracfarey


def plot_tune_diagram(
    orders: int | tuple[int] = (1, 2, 3),
    periodicity: int = 1,
    window: list = (0, 1, 0, 1),
    verbose: bool = False,
    legend: bool = False,
    show: bool = True,
    block: bool = False,
    debug: bool = False,
    axes: matplotlib.axes.Axes = None,
    linestyle: dict or str = None,
    linecolor: str or any = None,
    linewidth: int = None,
    addtolabel: str = None,
    **kwargs: dict[str, any],
) -> tuple[matplotlib.axes.Axes, list, list]:
    r"""
    Plot the tune diagram and resonance lines for the given *orders*, *periodicity*
    and *window*.

    The resonance equation is :math:`a\nu_x + b\nu_y = c`
    with :math:`a,b` and :math:`c` integers. The order is: :math:`N=abs(a)+abs(b)`.

    Parameters:
        orders: integer or tuple of integers larger than zero. Default (1, 2, 3).
        periodicity: periodicity of the machine, integer larger than zero. Default: 1.
        window: ``(nux_min, nux_max, nuy_min, nuy_max)``: tuple of 4 values for the
            tune minimum and maximum window. Default: (0, 1, 0, 1).
            *window* is ignored if the parameter axes is given.
        verbose: print verbose output.
        legend: print legend on the plot. Default: False.
        show: show plot. Default: True.
        block: passed to plot.show(). Default: False.
        debug: extra output to check line construction. Default: False.
        axes: :py:class:`~matplotlib.axes.Axes` for plotting the
            resonances. If :py:obj:`None`, a new figure will be created.
            Note that if *axes* are given then *window* is ignored.
        linestyle: line style for normal and skew resonances.

            If a dictionary is passed, it should contain
            {"normal": style1, "skew": style2}

            If 'dots' it uses dotted styles as linestyles. Equivalent to:
            {"normal": "dashdot", "skew": "dotted"}

            Default: {"normal": '-', "skew": '--'}
        linecolor: single color for all resonances. Default: custom palette.
            See :ref:`Lines Color and Width <color_width>`
        linewidth: line width for all resonances. Default: custom values.
            See :ref:`Lines Color and Width <color_width>`
        addtolabel: adds a string to the line label, e.g. for the fourth
            order normal resonance "4n"+addtolabel
    Keyword Args:
        only (str): if 'normal', plot only normal resonances.

            if 'skew' plots only skew resonances.

            Otherwise, ignored. See the notes on Normal and Skew convention.
        linedict (dict): dictionary of custom line styles. See notes below.

    Returns:
        Axes (matplotlib.axes.Axes): object from matplotlib.axes._axes
        legend_h (list):  list of handles for the legend
        legend_lab (list): list of labels for the legend

    NOTES:

    Normal and Skew convention:
    Line style is similar to reson.m from Matlab Middle Layer, MML, by L. Nadolski.
    Normal resonances are plotted with a continuous line.
    Skew resonances, i.e. N-abs(a) is odd, are plotted in dashed lines.

    .. _color_width:

    Lines Color and Width:
    Line style is similar to reson.m from Matlab Middle Layer, MML, by L. Nadolski.
    1st: black, width 4
    2nd: blue, width 3
    3rd: red, width 2
    4th: green, width 1
    5th: magenta, witdh 1
    6th: cyan, width 1
    7th: yellow, width 1
    8th: darkcyan, width 1
    9th: 'lightgreen', width 1
    10th to 15th: RGB increased in steps of 0.1, width 1

    Custom Style: You could pass a custom line style in a dictionary as
    ``linedict=mydictionary``, where mydictionary should contain two entries:
    dict("normal": normald, "skew": skewd).
    normald and skewd are also dictionaries, each entry contains as key
    the resonance order and as value the line properties to use in the plot.
    The default dictionary is created with :py:func:`create_linepalette`
    mydictionary = at.plot.resonances.create_linepalette()
    you could edit the needed entries.

    Raises:
        ValueError: if given resonances are lower than 0, or window is zero.
    """
    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    # print debugging output, equivalent to extra verbose
    debugprint = print if debug else lambda *a, **k: None

    # only plot normal or skew
    normalskew = [0, 1]
    onlyns = kwargs.pop("only", False)
    if onlyns == "normal":
        normalskew = [0]
    if onlyns == "skew":
        normalskew = [1]

    # orders could be a single int
    if isinstance(orders, int):
        orders = [orders]
    # check that all are larger than 0
    if sum(n < 1 for n in orders):
        raise ValueError("Negative resonances are not allowed.")

    if periodicity < 1:
        raise ValueError("Period must be a positive integer.")
    verboseprint(f"The periodicity is {periodicity}")

    verboseprint(f"The window is {window}")
    # check the window
    windowa = np.array(window)
    if windowa[0] == windowa[1]:
        raise ValueError("horizontal coordinates must be different")
    if windowa[2] == windowa[3]:
        raise ValueError("vertical coordinates must be different")
    if windowa[1] < windowa[0]:
        windowa[0], windowa[1] = windowa[1], windowa[0]
        warnings.warn("Swapping horizontal coordinates", stacklevel=1)
    if windowa[3] < windowa[2]:
        windowa[2], windowa[3] = windowa[3], windowa[2]
        warnings.warn("Swapping vertical coordinates", stacklevel=1)
    # get xlimits and ylimits
    the_axeslims = windowa.reshape((2, 2))

    maxreson2calc = np.max(orders)
    verboseprint(f"Farey max order={maxreson2calc}")

    # get the Farey collection, i.e., a list of farey sequences, one per order
    fareycollectionfloat = {}
    fareycollectionfrac = {}
    for nthorder in range(1, maxreson2calc + 1):
        farey, fracfarey = farey_sequence(nthorder)
        fareycollectionfloat[nthorder] = farey
        fareycollectionfrac[nthorder] = fracfarey
    verboseprint(f"the Farey collection is {fareycollectionfloat}")

    # plot configuration
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)
    else:
        verboseprint("Axes already exist, ignore window")
        the_axeslims = np.array([axes.get_xlim(), axes.get_ylim()])
    # min/max to plot lines with slopes
    minx = np.floor(the_axeslims[0, 0])
    minx = minx - periodicity - np.mod(minx, periodicity)
    maxx = np.ceil(the_axeslims[0, 1])
    maxx = maxx + periodicity - np.mod(maxx, periodicity)
    minmaxxdist = maxx - minx
    miny = np.floor(the_axeslims[1, 0])
    miny = miny - periodicity - np.mod(miny, periodicity)
    maxy = np.ceil(the_axeslims[1, 1])
    maxy = maxy + periodicity - np.mod(maxy, periodicity)
    minmaxydist = maxy - miny

    lprop = create_linepalette(
        linestyle=linestyle,
        linecolor=linecolor,
        linewidth=linewidth,
        addtolabel=addtolabel,
    )
    userprop = kwargs.pop("linedict", {})
    lprop["normal"].update(userprop.get("normal", {}))
    lprop["skew"].update(userprop.get("skew", {}))

    # we only need two points to define a line
    nauxpoints = 2

    # window min/max,horizontal and vertical
    axes.set_xlim(the_axeslims[0, :])
    axes.set_ylim(the_axeslims[1, :])
    # start to check the Farey collection, starting with 0
    collectaux1 = [0]
    idxtotype = {0: "normal", 1: "skew"}
    for nthorder in range(1, maxreson2calc + 1):
        debugprint(f"nthorder={nthorder}")
        collectaux2 = fareycollectionfloat[nthorder]
        thesteps = list(set(collectaux2) - set(collectaux1))
        debugprint(thesteps)
        chosenstep = min(thesteps)
        debugprint(f"chosenstep={chosenstep}")
        collectaux1 = collectaux2
        if nthorder in orders:
            # increase step by the periodicity in straight resonances
            debugprint("enter plotting horizontal straight lines")
            if 0 in normalskew:
                for iaux in np.arange(minx, maxx + 0.000001, periodicity * chosenstep):
                    axes.axvline(x=iaux, **lprop["normal"][nthorder])
            debugprint("enter plotting vertical straight lines")
            nsaux = np.mod(nthorder, 2)
            if nsaux in normalskew:
                for iaux in np.arange(miny, maxy + 0.000001, periodicity * chosenstep):
                    axes.axhline(y=iaux, **lprop[idxtotype[nsaux]][nthorder])
            # aeq*nux + beq*nuy = nthorder
            for aeq in range(1, nthorder):
                debugprint(f"enter plotting diagonals {aeq}")
                beq = nthorder - aeq
                chosenslope = 1.0 * aeq / beq
                ystep = 1.0 / beq
                debugprint(f"chosen slope={chosenslope}")
                debugprint(f"chosen ystep={ystep}")
                # calculate the variation on the vertical direction from window size
                yyaux = np.ceil(minmaxxdist * chosenslope + minmaxydist * ystep)
                y1aux = -yyaux + miny
                y2aux = yyaux + maxy
                # adapt to periodicity
                y1aux = periodicity * np.floor(y1aux / periodicity)
                y2aux = periodicity * np.ceil(y2aux / periodicity)
                debugprint(f"minx={minx},maxx={maxx},minmaxxdist={minmaxxdist}")
                debugprint(f"miny={miny},maxy={maxy},minmaxydist={minmaxydist}")
                debugprint(f"y1aux={y1aux},y2aux={y2aux}")
                xaux = np.linspace(minx, maxx, nauxpoints)
                debugprint(f"xaux={xaux}")
                nsaux = np.mod(beq, 2)
                if nsaux in normalskew:
                    for istep in np.arange(0, y2aux - y1aux + 0.0001, ystep):
                        y1line = (
                            -chosenslope * (xaux - minx) + y2aux - periodicity * istep
                        )
                        y2line = (
                            chosenslope * (xaux - minx) + y1aux + periodicity * istep
                        )
                        debugprint(f"y1line={y1line},y2line={y2line}")
                        axes.plot(xaux, y1line, **lprop[idxtotype[nsaux]][nthorder])
                        axes.plot(xaux, y2line, **lprop[idxtotype[nsaux]][nthorder])
    # include labels
    axes.set_xlabel(r"$\nu_x$")
    axes.set_ylabel(r"$\nu_y$")
    # printing legend if necessary
    handleall, labelall = axes.get_legend_handles_labels()
    hstyle = [hline._linestyle for hline in handleall]
    joinlabelstyle = [i + j for i, j in zip(labelall, hstyle)]
    uniquelabelstyle = sorted(set(joinlabelstyle))
    idxunique = [joinlabelstyle.index(uniqele) for uniqele in uniquelabelstyle]
    myleghandles = [handleall[idx] for idx in idxunique]
    myleglabels = [labelall[idx] for idx in idxunique]
    if legend:
        axes.legend(handles=myleghandles, labels=myleglabels, frameon=legend)
    if show:
        plt.show(block=block)

    return axes, myleghandles, myleglabels

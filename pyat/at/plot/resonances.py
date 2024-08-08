"""AT plotting functions related to resonances."""

import warnings
from fractions import Fraction

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy

__all__ = ["farey_sequence", "plot_tune_diagram"]

# 2024jul31 oblanco at ALBA CELLS


def create_linepalette(
    linestyle: str or dict = None,
    linecolor: str = None,
    linewidth: int = None,
) -> dict[str, any]:
    """
    Create a line palette to plot resonance lines.

    Arguments:
        linestyle: str or dictionary.
            If 'dots' it uses dotted styles as linestyles
                {"normal": "dashdot", "skew": "dotted"}
            If a dictionary is passed, it should contain
                {"normal": style1, "skew": style2}
            Default: {"normal": '-', "skew": '--'}
        linecolor: defines one color to be used. e.g. 'k'.
            Default: custom values. See :py:func:`plot_tune_diagram`
        linewidth: defines one integer value for the line width, e.g. 1.
            Default: custom values. See :py:func:`plot_tune_diagram`

    Returns:
        Dict dictionary contaning the line properties for resonance plots.
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
    # set the dictionary
    prop1n = {
        "color": mypalettecolor[1],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "1n",
    }
    prop2n = {
        "color": mypalettecolor[2],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[2],
        "label": "2n",
    }
    prop3n = {
        "color": mypalettecolor[3],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[3],
        "label": "3n",
    }
    prop4n = {
        "color": mypalettecolor[4],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[4],
        "label": "4n",
    }
    prop5n = {
        "color": mypalettecolor[5],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "5n",
    }
    prop6n = {
        "color": mypalettecolor[6],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "6n",
    }
    prop7n = {
        "color": mypalettecolor[7],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "7n",
    }
    prop8n = {
        "color": mypalettecolor[8],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "8n",
    }
    prop9n = {
        "color": mypalettecolor[9],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "9n",
    }
    prop10n = {
        "color": mypalettecolor[10],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "10n",
    }
    prop11n = {
        "color": mypalettecolor[11],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "11n",
    }
    prop12n = {
        "color": mypalettecolor[12],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "12n",
    }
    prop13n = {
        "color": mypalettecolor[13],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "13n",
    }
    prop14n = {
        "color": mypalettecolor[14],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "14n",
    }
    prop15n = {
        "color": mypalettecolor[15],
        "linestyle": mypalettestyle["normal"],
        "linewidth": mylinewidth[1],
        "label": "15n",
    }
    prop1s = {
        "color": mypalettecolor[1],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[1],
        "label": "1s",
    }
    prop2s = {
        "color": mypalettecolor[2],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[2],
        "label": "2s",
    }
    prop3s = {
        "color": mypalettecolor[3],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[3],
        "label": "3s",
    }
    prop4s = {
        "color": mypalettecolor[4],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[4],
        "label": "4s",
    }
    prop5s = {
        "color": mypalettecolor[5],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[5],
        "label": "5s",
    }
    prop6s = {
        "color": mypalettecolor[6],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[6],
        "label": "6s",
    }
    prop7s = {
        "color": mypalettecolor[7],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[7],
        "label": "7s",
    }
    prop8s = {
        "color": mypalettecolor[8],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[8],
        "label": "8s",
    }
    prop9s = {
        "color": mypalettecolor[9],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[9],
        "label": "9s",
    }
    prop10s = {
        "color": mypalettecolor[10],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[10],
        "label": "10s",
    }
    prop11s = {
        "color": mypalettecolor[11],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[11],
        "label": "11s",
    }
    prop12s = {
        "color": mypalettecolor[12],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[12],
        "label": "12s",
    }
    prop13s = {
        "color": mypalettecolor[13],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[14],
        "label": "13s",
    }
    prop14s = {
        "color": mypalettecolor[14],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[14],
        "label": "14s",
    }
    prop15s = {
        "color": mypalettecolor[15],
        "linestyle": mypalettestyle["skew"],
        "linewidth": mylinewidth[15],
        "label": "15s",
    }
    # assemble default dictionary by normal or skew
    propn = {
        1: prop1n,
        2: prop2n,
        3: prop3n,
        4: prop4n,
        5: prop5n,
        6: prop6n,
        7: prop7n,
        8: prop8n,
        9: prop9n,
        10: prop10n,
        11: prop11n,
        12: prop12n,
        13: prop13n,
        14: prop14n,
        15: prop15n,
    }
    props = {
        1: prop1s,
        2: prop2s,
        3: prop3s,
        4: prop4s,
        5: prop5s,
        6: prop6s,
        7: prop7s,
        8: prop8s,
        9: prop9s,
        10: prop10s,
        11: prop11s,
        12: prop12s,
        13: prop13s,
        14: prop14s,
        15: prop15s,
    }
    # assemble final dict
    return {"normal": propn, "skew": props}


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
        caux = numpy.floor((nthorder + bfarey) / dfarey) * cfarey - afarey
        daux = numpy.floor((nthorder + bfarey) / dfarey) * dfarey - bfarey
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
    orders: int or tuple = (1, 2, 3),
    periodicity: int = 1,
    window: list = (0, 1, 0, 1),
    verbose: bool = False,
    legend: bool = False,
    block: bool = False,
    debug: bool = False,
    axes: matplotlib.axes.Axes = None,
    **kwargs: dict[str, any],
) -> matplotlib.axes.Axes:
    r"""
    Plot the tune diagram and resonance lines for a given order, periodicity and window.

    Parameters:
        orders: integer or tuple of integers larger than zero. Default (1,2,3).
        periodicity: periodicity of the machine, integer larger than zero. Default: 1.
        window: (min_nux,max_nux,min_nuy,max_nuy) tuple of 4 values for the
            tune minimum and maximum window. Default: (0,1,0,1).
        verbose: print verbose output.
        legend: print legend on the plot. Default: False.
        block: passed to plot.show(). Default: False.
        debug: extra output to check line construction. Default: False.
        axes: :py:class:`~matplotlib.axes.Axes` for plotting the
            synoptic. If :py:obj:`None`, a new figure will be created. Otherwise,
            a new axes object sharing the same x-axis as the given one is created.
        kwargs:
            * only: if 'normal' plots only normal resonances.
                    if 'skew' plots only skew resonances.
                    Otherwise ignored.
            * linestyle: use it to pass a dictionary with custom line styles.
                See notes below.
            * linewidth = integer width to be used for all resonances.
                Default: See Color and Width below
            * linestyle: sets the line style for normal and skew resonances.
                If 'dots' is given it will use dashdot and dotted for normal
                    and skew resonances, respectively.
                A dictionary could be passed containing
                    mystyle = {"normal": style1, "skew":style2}
                    to plot using style1 and style2.
                Default: uses "-" and "--". See Normal and Skew convention.
            * linecolor: sets a single color for all the resonances.
                By default a custom palette is used. See Lines Color and Width.

    Returns:
        Axes object from matplotlib.axes._axes

    NOTES:
    The resonance equation is :math:`a\nu_x + b\nu_y = c`
    with :math:`a,b` and :math:`c` integers. The order :math:`N=abs(a)+abs(b)`.

    Normal and Skew convention:
    Line style is similar to reson.m from Matlab Middle Layer, MML, by L. Nadolski.
    Normal resonances are plotted with a continuous line.
    Skew resonances, i.e. N-abs(a) is odd, are plotted in dashed lines.

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

    Custom Style:
    You could pass a custom line style in a dictionary as
    linestyle = mydictionary,
    where mydictionary should contain two entries
    dict("normal": normald, "skew": skewd).
    normald and skewd are also dictionaries, each entry contains as key
    the resonance order and as value the line properties to use in the plot.

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
    windowa = numpy.array(window)
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

    maxreson2calc = numpy.max(orders)
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
        the_axeslims = numpy.array([axes.get_xlim(), axes.get_ylim()])
    # min/max to plot lines with slopes
    minx = numpy.floor(the_axeslims[0, 0])
    minx = minx - periodicity - numpy.mod(minx, periodicity)
    maxx = numpy.ceil(the_axeslims[0, 1])
    maxx = maxx + periodicity - numpy.mod(maxx, periodicity)
    minmaxxdist = maxx - minx
    miny = numpy.floor(the_axeslims[1, 0])
    miny = miny - periodicity - numpy.mod(miny, periodicity)
    maxy = numpy.ceil(the_axeslims[1, 1])
    maxy = maxy + periodicity - numpy.mod(maxy, periodicity)
    minmaxydist = maxy - miny

    linestyle = kwargs.pop("linestyle", None)
    linecolor = kwargs.pop("linecolor", None)
    linewidth = kwargs.pop("linewidth", None)
    defaultlprop = create_linepalette(
        linestyle=linestyle, linecolor=linecolor, linewidth=linewidth
    )
    lprop = kwargs.pop("linestyle", defaultlprop)

    # we only need to points to define a line
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
                for iaux in numpy.arange(
                    minx, maxx + 0.000001, periodicity * chosenstep
                ):
                    axes.axvline(x=iaux, **lprop["normal"][nthorder])
            debugprint("enter plotting vertical straight lines")
            nsaux = numpy.mod(nthorder, 2)
            if nsaux in normalskew:
                for iaux in numpy.arange(
                    miny, maxy + 0.000001, periodicity * chosenstep
                ):
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
                yyaux = numpy.ceil(minmaxxdist * chosenslope + minmaxydist * ystep)
                y1aux = -yyaux + miny
                y2aux = yyaux + maxy
                # adapt to periodicity
                y1aux = periodicity * numpy.floor(y1aux / periodicity)
                y2aux = periodicity * numpy.ceil(y2aux / periodicity)
                debugprint(f"minx={minx},maxx={maxx},minmaxxdist={minmaxxdist}")
                debugprint(f"miny={miny},maxy={maxy},minmaxydist={minmaxydist}")
                debugprint(f"y1aux={y1aux},y2aux={y2aux}")
                xaux = numpy.linspace(minx, maxx, nauxpoints)
                debugprint(f"xaux={xaux}")
                nsaux = numpy.mod(beq, 2)
                if nsaux in normalskew:
                    for istep in numpy.arange(0, y2aux - y1aux + 0.0001, ystep):
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
    myleghandles = []
    myleglabels = []
    handleall, labelall = axes.get_legend_handles_labels()
    myleglabels = sorted(set(labelall))
    idxunique = [labelall.index(mylab) for mylab in myleglabels]
    myleghandles = [handleall[idx] for idx in idxunique]
    if legend:
        axes.legend(handles=myleghandles, labels=myleglabels, frameon=legend)
    plt.show(block=block)

    return axes

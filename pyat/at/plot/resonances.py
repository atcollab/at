"""AT plotting functions related to resonances."""

import warnings
from fractions import Fraction

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy
from matplotlib.figure import Figure

__all__ = ["farey_sequence", "plot_tune2d_resonances"]

# 2024jul31 oblanco at ALBA CELLS

# create default dictionary with line properties
widthmod = 5  # should it be a variable??? it limits the linewidth
mypalettecolor = {
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
mypalettestyle = {"normal": "-", "skew": "--"}
prop1n = {
    "color": mypalettecolor[1],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(4, widthmod),
    "label": "1n",
}
prop2n = {
    "color": mypalettecolor[2],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(3, widthmod),
    "label": "2n",
}
prop3n = {
    "color": mypalettecolor[3],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(2, widthmod),
    "label": "3n",
}
prop4n = {
    "color": mypalettecolor[4],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(2, widthmod),
    "label": "4n",
}
prop5n = {
    "color": mypalettecolor[5],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "5n",
}
prop6n = {
    "color": mypalettecolor[6],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "6n",
}
prop7n = {
    "color": mypalettecolor[7],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "7n",
}
prop8n = {
    "color": mypalettecolor[8],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "8n",
}
prop9n = {
    "color": mypalettecolor[9],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "9n",
}
prop10n = {
    "color": mypalettecolor[10],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "10n",
}
prop11n = {
    "color": mypalettecolor[11],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "11n",
}
prop12n = {
    "color": mypalettecolor[12],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "12n",
}
prop13n = {
    "color": mypalettecolor[13],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "13n",
}
prop14n = {
    "color": mypalettecolor[14],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "14n",
}
prop15n = {
    "color": mypalettecolor[15],
    "linestyle": mypalettestyle["normal"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "15n",
}
prop1s = {
    "color": mypalettecolor[1],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(4, widthmod),
    "label": "1s",
}
prop2s = {
    "color": mypalettecolor[2],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(3, widthmod),
    "label": "2s",
}
prop3s = {
    "color": mypalettecolor[3],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(2, widthmod),
    "label": "3s",
}
prop4s = {
    "color": mypalettecolor[4],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(2, widthmod),
    "label": "4s",
}
prop5s = {
    "color": mypalettecolor[5],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "5s",
}
prop6s = {
    "color": mypalettecolor[6],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "6s",
}
prop7s = {
    "color": mypalettecolor[7],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "7s",
}
prop8s = {
    "color": mypalettecolor[8],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "8s",
}
prop9s = {
    "color": mypalettecolor[9],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "9s",
}
prop10s = {
    "color": mypalettecolor[10],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "10s",
}
prop11s = {
    "color": mypalettecolor[11],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "11s",
}
prop12s = {
    "color": mypalettecolor[12],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "12s",
}
prop13s = {
    "color": mypalettecolor[13],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "13s",
}
prop14s = {
    "color": mypalettecolor[14],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
    "label": "14s",
}
prop15s = {
    "color": mypalettecolor[15],
    "linestyle": mypalettestyle["skew"],
    "linewidth": numpy.mod(1, widthmod),
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
defaultlprop = {0: propn, 1: props}


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


def plot_tune2d_resonances(
    orders: int or list = [1, 2, 3],
    period: int = 1,
    window: list = [0, 1, 0, 1],
    verbose: bool = False,
    legend: bool = False,
    block: bool = False,
    debug: bool = False,
    **kwargs: dict[str, any],
) -> Figure:
    r"""
    Plot the tune 2D resonances for a given order, period and window.

    Parameters:
        orders: integer or list of integers larger than zero. Default: [1,2,3]
        period: integer larger than zero; periodicity of the machine. Default: 1.
        window: [min_nux,max_nux,min_nuy,max_nuy] list of 4 values for the
                tune minimum and maximum window. Default:[0,1,0,1].
        verbose: print verbose output.
        legend: print legend on the plot. Default: False.
        block: passed to plot.show(). Default: False.
        debug: extra output to check line construction. Default: False.
        kwargs:
            * onlyns: if 'n' plots only normal resonances.
                      if 's' plots only skew resonances.
                      Otherwise ignored.
            * linestyle: use it to pass a dictionary with custom line styles.
                See notes below.

    Returns:
        Figure object

    NOTES:
    The resonance equation is :math:`a\nu_x + b\nu_y = c`
    with :math:`a,b` and :math:`c` integers. The order :math:`N=abs(a)+abs(b)`.

    Normal and Skew convention:
    Line style is similar to reson.m from Matlab Middle Layer, MML, by L. Nadolski.
    Normal resonances are plotted with a continuous line.
    Skew resonances, i.e. N-abs(a) is odd, are plotted in dashed lines.

    Color and thickness:
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
    dict(0: normal, 1: skew).
    normal and skew are also dictionaries, each entry contains the
    line properties of the nth resonance to plot.

    Raises:
        ValueError: if given resonances are lower than 0, or window is zero.
    """
    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    # print debugging output, equivalent to extra verbose
    debugprint = print if debug else lambda *a, **k: None

    # only plot normal or skew
    normalskew = [0, 1]
    onlyns = kwargs.pop("onlyns", False)
    if onlyns == "n":
        normalskew = [0]
    if onlyns == "s":
        normalskew = [1]

    # orders could be a single int
    if isinstance(orders, int):
        orders = [orders]
    # check that all are larger than 0
    if sum(n < 1 for n in orders):
        raise ValueError("Negative resonances are not allowed.")

    if period < 1:
        raise ValueError("Period must be a positive integer.")
    verboseprint(f"The period is {period}")

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

    # min/max to plot lines with slopes
    minx = numpy.floor(the_axeslims[0, 0])
    minx = minx - period - numpy.mod(minx, period)
    maxx = numpy.ceil(the_axeslims[0, 1])
    maxx = maxx + period - numpy.mod(maxx, period)
    minmaxxdist = maxx - minx
    miny = numpy.floor(the_axeslims[1, 0])
    miny = miny - period - numpy.mod(miny, period)
    maxy = numpy.ceil(the_axeslims[1, 1])
    maxy = maxy + period - numpy.mod(maxy, period)
    minmaxydist = maxy - miny

    lprop = kwargs.pop("linestyle", defaultlprop)

    # we only need to points to define a line
    nauxpoints = 2

    # plot configuration
    fig = plt.figure()
    # window min/max,horizontal and vertical
    plt.xlim(the_axeslims[0, :])
    plt.ylim(the_axeslims[1, :])
    # start to check the Farey collection, starting with 0
    collectaux1 = [0]
    for nthorder in range(1, maxreson2calc + 1):
        debugprint(f"nthorder={nthorder}")
        collectaux2 = fareycollectionfloat[nthorder]
        thesteps = list(set(collectaux2) - set(collectaux1))
        debugprint(thesteps)
        chosenstep = min(thesteps)
        debugprint(f"chosenstep={chosenstep}")
        collectaux1 = collectaux2
        if nthorder in orders:
            # increase step by the period in straight resonances
            debugprint("enter plotting horizontal straight lines")
            if 0 in normalskew:
                for iaux in numpy.arange(minx, maxx + 0.000001, period * chosenstep):
                    plt.axvline(x=iaux, **lprop[0][nthorder])
            debugprint("enter plotting vertical straight lines")
            nsaux = numpy.mod(nthorder, 2)
            if nsaux in normalskew:
                for iaux in numpy.arange(miny, maxy + 0.000001, period * chosenstep):
                    plt.axhline(y=iaux, **lprop[nsaux][nthorder])
            for iaux in range(1, nthorder):
                debugprint(f"enter plotting diagonals {iaux}")
                aeq = iaux
                beq = nthorder - aeq
                chosenslope = 1.0*aeq / beq
                diagstep = 1.0 / beq
                debugprint(f"chosen slope={chosenslope}")
                debugprint(f"chosen diagstep={diagstep}")
                # get the vertical limits with diagonals
                a1aux = (
                    -numpy.ceil(minmaxxdist*chosenslope + minmaxydist*diagstep)  + miny
                )
                a2aux = (
                    numpy.ceil(minmaxxdist*chosenslope + minmaxydist*diagstep) + maxy
                )
                debugprint(f"minx={minx},maxx={maxx},minmaxxdist={minmaxxdist}")
                debugprint(f"miny={miny},maxy={maxy},minmaxydist={minmaxydist}")
                debugprint(f"a1aux={a1aux},a2aux={a2aux}")
                xaux = numpy.linspace(minx, maxx, nauxpoints)
                debugprint(f"xaux={xaux}")
                nsaux = numpy.mod(beq, 2)
                if nsaux in normalskew:
                    for istep in numpy.arange(0, a2aux - a1aux + 0.0001, diagstep):
                        y1line = -chosenslope * (xaux - minx) + a2aux - period * istep
                        y2line =  chosenslope * (xaux - minx) + a1aux + period * istep
                        debugprint(f"y1line={y1line},y2line={y2line}")
                        plt.plot(xaux, y1line, **lprop[nsaux][nthorder])
                        plt.plot(xaux, y2line, **lprop[nsaux][nthorder])
    # include labels
    plt.xlabel(r"$\nu_x$")
    plt.ylabel(r"$\nu_y$")
    # printing legend if necessary
    myleghandles = []
    if legend:
        for nthorder in orders:
            for nsaux in normalskew:
                dictaux = lprop[nsaux][nthorder]
                myleghandles.append(mlines.Line2D([], [], **dictaux))
    plt.legend(handles=myleghandles, frameon=legend)
    plt.show(block=block)

    return fig

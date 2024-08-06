"""AT plotting functions related to resonances."""

from fractions import Fraction
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy

from ..lattice.utils import AtError, AtWarning

__all__ = ["farey_sequence", "plot_tune2D_resonances"]


def farey_sequence(nthorder: int, verbose: bool = False):
    """
    Returns the Farey sequence, and the resonance sequence of nth order.

    Arguments:
        nthorder: natural number bigger than 0
    Options:
        verbose: prints extra info. Default: False
    Returns:
        fareyseqfloat: list of elements with the Farey sequence in Float format.
            See Eqs.(1,2,3) of [1].
        fareyseqfrac: list of elements with the Farey sequence in Fraction format.
            See Eqs.(1,2,3) of [1].

    [1] R.Tomas. 'From Farey sequences to resonance diagrams. Phys.Rev.Acc.Beams 17, 014001 (2014)'
    """
    verboseprint = print if verbose else lambda *a, **k: None
    verboseprint(f"nthorder={nthorder}")

    farey = []
    fracfarey = []
    af = 0
    bf = 1
    cf = 1
    df = nthorder
    farey.append(0)
    farey.append(1 / df)
    fracfarey.append(Fraction(0))
    fracfarey.append(Fraction(1, df))
    idx = 0
    while (farey[-1] < 1) and (idx < 100):
        idx += 1
        caux = numpy.floor((nthorder + bf) / df) * cf - af
        daux = numpy.floor((nthorder + bf) / df) * df - bf
        af = cf
        bf = df
        cf = int(caux)
        df = int(daux)
        farey.append(cf / df)
        fracfarey.append(Fraction(cf, df))
    verboseprint(f"farey_float{nthorder}= {farey}")
    verboseprint(f"farey_frac{nthorder} = {fracfarey}")
    return farey, fracfarey


def plot_tune2D_resonances(
    orders: any = [1, 2, 3],
    period: int = 1,
    window: list = [0, 1, 0, 1],
    verbose: bool = False,
    **kwargs: dict[any],
):
    r"""
    This function plots the tune 2D resonances for a given order, period and window.

    Options:
    orders: integer or list of integers larger than zero. Default: [1,2,3]
    period: integer larger than zero; periodicity of the machine. Default: 1
    window: [min_nux,max_nux,min_nuy,max_nuy] list of 4 values for the tune minimum
      and maximum window. Default:[0,1,0,1]
    verbose: print verbose output
    includelegend: print legend on the plot. Default: False
    onlyns: if 'n' plots only normal resonances.
            if 's' plots only skew resonances.
            Otherwise ignored.
    custom_linesty: use it to pass a dictionary with custom line styles. See notes
      below.

    NOTES:
    The resonance equation is:
      :math:`a\nu_x + b\nu_y = c`
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
      custum_linesty = mydictionary
    where mydictionary should contain two entries, dict(0: prop_normal, 1: prop_skew)
      where prop_normal and prop_skew are also dictionaries starting at zero,
      i.e. the_nth_order-1. Each entry contains the line properties of the nth-1 resonance to plot.
    """
    # 2024jul31 oblanco at ALBA CELLS

    # verboseprint to check flag only once
    verboseprint = print if verbose else lambda *a, **k: None

    # print debugging output, equivalent to extra verbose
    debugverbose = kwargs.pop("debugverbose", False)
    debugprint = print if debugverbose else lambda *a, **k: None

    # block the standard output in terminal when plotting
    block = kwargs.pop("block", False)

    # plot with legend
    includelegend = kwargs.pop("includelegend", False)

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
    listresonancestoplot = list(numpy.array(orders) - 1)
    # check that all are larger than 0
    if sum(n < 0 for n in listresonancestoplot):
        AtError("Negative resonances are not allowed")

    theperiod = period
    verboseprint(f"The period is {theperiod}")
    verboseprint(f"The window is {window}")

    # check the window
    windowa = numpy.array(window)
    if windowa[0] == windowa[1]:
        AtError("horizontal coordinates must be different")
    if windowa[2] == windowa[3]:
        AtError("vertical coordinates must be different")
    if windowa[1] < windowa[0]:
        AtWarning("Swapping horizontal coordinates")
        windowa[0], windowa[1] = windowa[1], windowa[0]
    if windowa[3] < windowa[2]:
        AtWarning("Swapping vertical coordinates")
        windowa[2], windowa[3] = windowa[3], windowa[2]
    # get xlimits and ylimits
    the_axeslims = windowa.reshape((2, 2))

    # horizontal and vertical borders
    borders = numpy.eye(2)

    maxreson2calc = numpy.max(listresonancestoplot) + 1
    verboseprint(f"Farey max order={maxreson2calc}")

    # get the Farey collection, i.e., a list of farey sequences, one per order
    fareycollectionfloat = []
    fareycollectionfrac = []
    for nthorder in range(1, maxreson2calc + 1):
        farey, fracfarey = farey_sequence(nthorder)
        fareycollectionfloat.append(farey.copy())
        fareycollectionfrac.append(fracfarey.copy())
    verboseprint(f"the Farey collection is {fareycollectionfloat}")

    # plot configuration
    fig = plt.figure()
    # window min/max,horizontal and vertical
    plt.xlim(the_axeslims[0, :])
    plt.ylim(the_axeslims[1, :])
    # min/max to plot lines with slopes
    minx = numpy.floor(the_axeslims[0, 0])
    minx = minx - theperiod - numpy.mod(minx, theperiod)
    maxx = numpy.ceil(the_axeslims[0, 1])
    maxx = maxx + theperiod - numpy.mod(maxx, theperiod)
    minmaxxdist = maxx - minx
    miny = numpy.floor(the_axeslims[1, 0])
    miny = miny - theperiod - numpy.mod(miny, theperiod)
    maxy = numpy.ceil(the_axeslims[1, 1])
    maxy = maxy + theperiod - numpy.mod(maxy, theperiod)
    minmaxydist = maxy - miny
    minmaxdist = max([minmaxxdist, minmaxydist])

    # dictionary with line properties
    widthmod = 5  # ??? maybe a variable
    mypalettecolor = {
        0: "k",
        1: "b",
        2: "r",
        3: "g",
        4: "m",
        5: "c",
        6: "y",
        7: "darkcyan",
        8: "lightgreen",
        9: (0.1, 0.1, 0.1),
        10: (0.1, 0.1, 0.1),
        11: (0.2, 0.2, 0.2),
        12: (0.3, 0.3, 0.3),
        13: (0.4, 0.4, 0.4),
        14: (0.5, 0.5, 0.5),
        15: (0.6, 0.6, 0.6),
    }
    mypalettestyle = {0: "-", 1: "--"}
    prop1n = dict(
        color=mypalettecolor[0],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(4, widthmod),
        label="1n",
    )
    prop2n = dict(
        color=mypalettecolor[1],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(3, widthmod),
        label="2n",
    )
    prop3n = dict(
        color=mypalettecolor[2],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(2, widthmod),
        label="3n",
    )
    prop4n = dict(
        color=mypalettecolor[3],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(2, widthmod),
        label="4n",
    )
    prop5n = dict(
        color=mypalettecolor[4],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="5n",
    )
    prop6n = dict(
        color=mypalettecolor[5],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="6n",
    )
    prop7n = dict(
        color=mypalettecolor[6],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="7n",
    )
    prop8n = dict(
        color=mypalettecolor[7],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="8n",
    )
    prop9n = dict(
        color=mypalettecolor[8],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="9n",
    )
    prop10n = dict(
        color=mypalettecolor[9],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="10n",
    )
    prop11n = dict(
        color=mypalettecolor[10],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="11n",
    )
    prop12n = dict(
        color=mypalettecolor[11],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="12n",
    )
    prop13n = dict(
        color=mypalettecolor[12],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="13n",
    )
    prop14n = dict(
        color=mypalettecolor[13],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="14n",
    )
    prop15n = dict(
        color=mypalettecolor[14],
        linestyle=mypalettestyle[0],
        linewidth=numpy.mod(1, widthmod),
        label="15n",
    )
    prop1s = dict(
        color=mypalettecolor[0],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(4, widthmod),
        label="1s",
    )
    prop2s = dict(
        color=mypalettecolor[1],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(3, widthmod),
        label="2s",
    )
    prop3s = dict(
        color=mypalettecolor[2],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(2, widthmod),
        label="3s",
    )
    prop4s = dict(
        color=mypalettecolor[3],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(2, widthmod),
        label="4s",
    )
    prop5s = dict(
        color=mypalettecolor[4],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="5s",
    )
    prop6s = dict(
        color=mypalettecolor[5],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="6s",
    )
    prop7s = dict(
        color=mypalettecolor[6],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="7s",
    )
    prop8s = dict(
        color=mypalettecolor[7],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="8s",
    )
    prop9s = dict(
        color=mypalettecolor[8],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="9s",
    )
    prop10s = dict(
        color=mypalettecolor[9],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="10s",
    )
    prop11s = dict(
        color=mypalettecolor[10],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="11s",
    )
    prop12s = dict(
        color=mypalettecolor[11],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="12s",
    )
    prop13s = dict(
        color=mypalettecolor[12],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="13s",
    )
    prop14s = dict(
        color=mypalettecolor[13],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="14s",
    )
    prop15s = dict(
        color=mypalettecolor[14],
        linestyle=mypalettestyle[1],
        linewidth=numpy.mod(1, widthmod),
        label="15s",
    )
    # assemble default dictionary by normal or skew
    propn = {
        0: prop1n.copy(),
        1: prop2n.copy(),
        2: prop3n.copy(),
        3: prop4n.copy(),
        4: prop5n.copy(),
        5: prop6n.copy(),
        6: prop7n.copy(),
        7: prop8n.copy(),
        8: prop9n.copy(),
        9: prop10n.copy(),
        10: prop11n.copy(),
        11: prop12n.copy(),
        12: prop13n.copy(),
        13: prop14n.copy(),
        14: prop15n.copy(),
    }
    props = {
        0: prop1s.copy(),
        1: prop2s.copy(),
        2: prop3s.copy(),
        3: prop4s.copy(),
        4: prop5s.copy(),
        5: prop6s.copy(),
        6: prop7s.copy(),
        7: prop8s.copy(),
        8: prop9s.copy(),
        9: prop10s.copy(),
        10: prop11s.copy(),
        11: prop12s.copy(),
        12: prop13s.copy(),
        13: prop14s.copy(),
        14: prop15s.copy(),
    }
    # assemble final dict
    defaultlprop = {0: propn.copy(), 1: props.copy()}
    lprop = kwargs.pop("customlinesty", defaultlprop)

    # we only need to points to define a line
    nauxpoints = 2

    # start to check the Farey collection, starting with 0
    collectaux1 = [0]
    for nthorder in range(0, maxreson2calc):
        debugprint(f"nthorder={nthorder}")
        collectaux2 = fareycollectionfloat[nthorder]
        thesteps = list(set(collectaux2) - set(collectaux1))
        debugprint(thesteps)
        chosenstep = min(thesteps)
        debugprint(f"chosenstep={chosenstep}")
        collectaux1 = collectaux2
        if not (nthorder in listresonancestoplot):
            continue
        else:
            # increase step by the period in straight resonances
            debugprint("enter plotting horizontal straight lines")
            if 0 in normalskew:
                for iaux in numpy.arange(minx, maxx + 0.000001, theperiod * chosenstep):
                    plt.axvline(x=iaux, **lprop[0][nthorder].copy())
            debugprint("enter plotting vertical straight lines")
            nsaux = numpy.mod(nthorder + 1, 2)
            if nsaux in normalskew:
                for iaux in numpy.arange(miny, maxy + 0.000001, theperiod * chosenstep):
                    plt.axhline(y=iaux, **lprop[nsaux][nthorder].copy())
            for iaux in range(nthorder):
                debugprint(f"enter plotting diagonals {iaux}")
                aeq = iaux + 1
                beq = nthorder + 1 - aeq
                chosenslope = -aeq / beq
                diagstep = 1 / beq
                debugprint(f"chosen slope={chosenslope}")
                debugprint(f"chosen diagstep={diagstep}")
                # get the vertical limits with diagonals
                a1 = (
                    -numpy.ceil(2 * (minmaxxdist + minmaxydist) / (-chosenslope)) + miny
                )
                a2 = numpy.ceil(2 * (minmaxxdist + minmaxydist) / (-chosenslope)) + maxy
                debugprint(f"minx={minx},maxx={maxx},minmaxxdist={minmaxxdist}")
                debugprint(f"miny={miny},maxy={maxy},minmaxydist={minmaxydist}")
                debugprint(f"a1={a1},a2={a2}")
                xaux = numpy.linspace(minx, maxx, nauxpoints)
                debugprint(f"xaux={xaux}")
                nsaux = numpy.mod(beq, 2)
                if nsaux in normalskew:
                    for istep in numpy.arange(0, a2 - a1 + 0.0001, diagstep):
                        y1line = chosenslope * (xaux - minx) + theperiod * istep + miny
                        y2line = -chosenslope * (xaux - minx) - theperiod * istep + maxy
                        debugprint(f"y1line={y1line},y2line={y2line}")
                        plt.plot(xaux, y1line, **lprop[nsaux][nthorder])
                        plt.plot(xaux, y2line, **lprop[nsaux][nthorder])
    # include labels
    plt.xlabel(r"$\nu_x$")
    plt.ylabel(r"$\nu_y$")
    # printing legend if necessary
    myleghandles = []
    if includelegend:
        for nthorder in listresonancestoplot:
            for nsaux in normalskew:
                dictaux = lprop[nsaux][nthorder].copy()
                myleghandles.append(mlines.Line2D([], [], **dictaux))
    plt.legend(handles=myleghandles, frameon=includelegend)
    plt.show(block=block)

    return fig

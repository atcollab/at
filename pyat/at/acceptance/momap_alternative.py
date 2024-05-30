import numpy
import time

_all_ = ["momaperture_project2start","projectrefpts"]

# 2024may29 oblanco at ALBA CELLS. First working version
#           Based on the MATLAB implementation by Z.Marti at ALBA


def momaperture_project2start(ring, *args, **kwargs):
    """
    :py:func:`momap_project2start` calculates the local momemtum aperture.
    It is a binary search of the negative and positive energy thresholds
    of stability around the closed orbit.

    For a given energy offset the particles are first tracked from every
    reference point to the end of the ring, and then, only particles with
    differing 6D coordinates are tracked together.  The surviving particles
    continue the boundary search with a new energy until the boundary is
    found or the limit of the rf bucket is reached.

    Usage:
      >>> momaperture_project2start(ring)

    Parameters:
        ring: list of elements

    Keyword arguments:
      refpts: Selects the location of coordinates output.
        See ":ref:`Selecting elements in a lattice <refpts>`"
      nturns: number of turns to be tracked. Default 1000
      dptol: energy offset resolution. Default 1e-4
      add_offset: (2, N) offsets to be added to the transverse coordinates
        on the N reference points. Default 1e-5 m
      eu_ini: maximum energy boundary. Default: the rf bucket height
      orbit: (N,6) offsets to be added on the N reference points.
        Default, the closed orbit
      verbose: print in the standard output additional info. Default False
      use_mp (bool): Default True. See :py:func:.`lattice_track`
      omp_num_threads (int): See :py:fun:`lattice_track`
      pool_size: See :py:func:`lattice_track`
      start_method: See :py:fun:`lattice_track`

    Returns:
      dnp: (2,N) array with negative and positive stable energy boundaries
        for the N reference points

    ..note::
      * This function does a quick search, but, it is succeptible to miss
        islands of instability due to the varying energy step.
    """
    # verboseprint to check flag only once
    verbose = False
    if "verbose" in kwargs:
        verbose = bool(kwargs["verbose"])
    verboseprint = print if verbose else lambda *a, **k: None

    if "refpts" in kwargs:
        rp = kwargs["refpts"]
    else:
        rp = ring.uint32_refpts(range(len(ring)))
    nrp = len(rp)
    verboseprint(f"Using {nrp} reference points")

    if "nturns" in kwargs:
        nturns = int(kwargs["nturns"])
    else:
        nturns = 1000
    verboseprint(f"Track over {nturns} turns")

    if "dptol" in kwargs:
        dptol = float(kwargs["dptol"])
    else:
        dptol = 1e-4
    verboseprint(f"Energy resolution {dptol}")

    if "add_offset" in kwargs:
        add_offset = kwargs["add_offset"]
        verboseprint(f"Add user offsets")
    else:
        dxy = 1e-5
        add_offset = numpy.tile(dxy, [2, nrp])
        verboseprint(f"Adding default transverse offsets {dxy}")

    # get the tracking options
    dicttrack = {}
    lpartrack = ["omp_num_threads", "use_mp", "pool_size", "start_method"]
    for pp in lpartrack:
        if pp in kwargs:
            dicttrack.update({pp: kwargs[pp]})
    # default to parallel
    if not "use_mp" in dicttrack:
        dicttrack.update({"use_mp": True})

    # use radiation parameters to get the rf bucket
    pars = ring.radiation_parameters()

    # get energy bucket (db) max height.
    # S.Y.Lee, 4th Edition, Eqs 3.37, 3.38, Sect. II.2 Bucket area
    h = ring.harmonic_number
    etac = pars.etac
    E = pars.E0
    phis = pars.phi_s
    betar = ring.beta
    V = ring.get_rf_voltage()
    Y = numpy.sqrt(
        numpy.abs(numpy.cos(phis) - 0.5 * (numpy.pi - 2 * phis) * numpy.sin(phis))
    )
    db = numpy.sqrt(2 * V / (numpy.pi * betar**2 * E * h * numpy.abs(etac))) * Y
    verboseprint(f"Bucket height {db}")

    # first guess
    es_ini = 0
    et_ini = db / 2
    if "eu_ini" in kwargs:
        eu_ini = kwargs["eu_ini"]
        verboseprint(f"Using the users max boundary {eu_ini}")
    else:
        eu_ini = db
        verboseprint("Using the bucket height as maximum boundary")

    if "orbit" in kwargs:
        orbit_s = kwargs["orbit"]
        verboseprint("Using the users orbit")
    else:
        _, orbit_s = ring.find_orbit(rp)
        verboseprint("Using the closed orbit")

    # start scan
    # eu, unstable energy
    # es, stable energy
    # et, test energy
    etnp = numpy.ones((2, nrp))
    thesign = [-1, 1]  # first negative side, then positive side
    verbosesign = ["negative", "positive"]
    for i in [0, 1]:
        verboseprint(f"Search on the {verbosesign[i]} side")
        et = thesign[i] * et_ini * numpy.ones(nrp)
        eu = thesign[i] * eu_ini * numpy.ones(nrp)
        es = thesign[i] * es_ini * numpy.ones(nrp)
        de = 1
        iteration = 0
        while iteration < 100 and (de > dptol):  # safety limit on iterations
            t0 = time.time()
            # L is a mask, True for lost particles
            L = multirefpts_track_islost(
                ring, rp, et, orbit_s, nturns, add_offset, **dicttrack
            )
            es[~L] = et[~L]
            eu[L] = et[L]
            et = (es + eu) / 2
            de = max(abs(es - eu))
            iteration = iteration + 1
            outmsg = (
                f"Iteration {iteration} in {verbosesign[i]} side",
                f" took {format(time.time()-t0):.3} s.",
                f" {de=}, {dptol=}",
            )
            verboseprint("".join(outmsg))
        etnp[i, :] = et
    etnp = etnp.T

    return etnp

def projectrefpts(ring, startrefpts, **kwargs):
    """
    This function tracks the particles from multiple refpts to a single
    refpts.
    """
    lenring = len(ring)
    srps = startrefpts
    nrps = len(srps)

    # verboseprint to check flag only once
    erp = lenring
    verbose = False
    if "verbose" in kwargs:
        verbose = bool(kwargs["verbose"])
    verboseprint = print if verbose else lambda *a, **k: None

    if 'endrefpt' in kwargs:
        erp = kwargs['endrefpt']
        verboseprint(f'Project particles to start of element index {erp}')
    else:
        erp = lenring
        verboseprint('Project to end point')

    if "particles" in kwargs:
        orbit_s = kwargs["particles"]
        verboseprint("Using the users particles")
    else:
        _, orbit_s = ring.find_orbit(rp)
        verboseprint("Using the closed orbit")
        dxy = 1e-5
        add_offset = numpy.tile(dxy, [2, nrp])
        verboseprint(f"Adding default transverse offsets {dxy}")

    zin = numpy.zeros((6, nrps))
    zout = numpy.zeros((6, nrps))
    lostpart = numpy.ones((nrps), dtype=bool)
    # first, track the remaining portion of the ring
    for i in range(nrps):
        ring_downstream = ring.rotate(srps[i])
        zin[:, i] = orbit_s[i, :].copy()
        zo, _, dout1 = ring_downstream.track(
            zin[:, i], nturns=1, refpts=lenring - rp[i], losses=True
        )
        lostpart[i] = dout1["loss_map"]["islost"][0]
        zout[:, i] = zo[:, 0, 0, 0]
    return zout, lostpart


def multirefpts_track_islost(ring, refpts, e, orbit, nturns, initcoord, **dicttrack):
    """
    Returns a boolean array: tells whether the particle launched is lost.
    True means lost
    """

    lenring = len(ring)
    rp = refpts
    nrp = len(rp)
    lostpart = numpy.ones((nrp), dtype=bool)
    zin = numpy.zeros((6, nrp))
    zout = numpy.zeros((6, nrp))
    issmall = 1e-6
    eps = numpy.finfo(float).eps
    istiny = 100 * eps

    # first, track the remaining portion of the ring
    for i in range(nrp):
        ring_downstream = ring.rotate(rp[i])  # [rp[i]:]
        zin[:, i] = orbit[i, :].copy()
        zin[0, :] = zin[0, :] + initcoord[0, :]
        zin[2, :] = zin[2, :] + initcoord[1, :]
        zin[4, :] = zin[4, :] + e
        zo, _, dout1 = ring_downstream.track(
            zin[:, i], nturns=1, refpts=lenring - rp[i], losses=True
        )
        lostpart[i] = dout1["loss_map"]["islost"][0]
        zout[:, i] = zo[:, 0, 0, 0]

    cntalive = nrp - sum(lostpart)
    aliveatringend = ~lostpart
    zinaliveaux = numpy.squeeze(zout[:, aliveatringend])
    zinalive = numpy.asfortranarray(zinaliveaux.copy())

    # second, track particles that have survived the ring to the end
    if cntalive == 1:
        zaux1, _, dout2 = ring.track(
            zinalive, nturns=nturns, refpts=len(ring), losses=True
        )
        lostpart[aliveatringend] = dout2["loss_map"]["islost"][0]
    elif cntalive > 1:
        # search for non numerically similar (istiny = 100 times eps) particles
        closenessmatrix = numpy.zeros((cntalive, cntalive), dtype=bool)
        for i in range(cntalive):
            closenessmatrix[i, :] = [
                numpy.allclose(zinalive[:, j], zinalive[:, i], atol=istiny)
                for j in range(cntalive)
            ]
        _, rowidx = numpy.indices((cntalive, cntalive))
        maxidx = numpy.max(rowidx * closenessmatrix, 1)
        uniqueidx, _, rep_index = numpy.unique(
            maxidx, return_index=True, return_inverse=True
        )
        # dummy track to later reuse the ring
        zaux2 = numpy.array([initcoord[0, 0], 0, initcoord[1, 0], 0, issmall, 0]).T
        ring.track(zaux2, nturns=1)
        # track non-numerically similar particles
        zaux3, _, dout3 = ring.track(
            zinalive[:, uniqueidx],
            nturns=nturns,
            keep_lattice=True,
            losses=True,
            **dicttrack,
        )
        replossmask = dout3["loss_map"]["islost"]
        # copy result for numerically similar (100 times eps) particles
        losteaux = replossmask[rep_index]
        # now, group the particles that did not pass the first turn
        lostpart[aliveatringend] = losteaux

    return lostpart

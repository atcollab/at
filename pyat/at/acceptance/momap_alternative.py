"""momap_alternative."""

import time

import numpy

__all__ = ["momaperture_project2start", "projectrefpts"]

# 2024may29 oblanco at ALBA CELLS. First working version
#           Based on the MATLAB implementation by Z.Marti at ALBA


def momaperture_project2start(ring: list, **kwargs: dict[str, any]) -> numpy.ndarray:
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

    rps = kwargs.pop("refpts", ring.uint32_refpts(range(len(ring))))
    nrps = len(rps)
    verboseprint(f"Using {nrps} reference points")

    nturns = kwargs.pop("nturns", 1000)
    verboseprint(f"Track over {nturns} turns")

    dptol = kwargs.pop("dptol", 1e-4)
    verboseprint(f"Energy resolution {dptol}")

    if "add_offset" in kwargs:
        add_offset = kwargs["add_offset"]
        verboseprint("Add user offsets")
    else:
        dxy = 1e-5
        add_offset = numpy.tile(dxy, [2, nrps])
        verboseprint(f"Adding default transverse offsets {dxy}")

    # get the tracking options
    dicttrack = {}
    lpartrack = ["omp_num_threads", "use_mp", "pool_size", "start_method"]
    for theparameter in lpartrack:
        if theparameter in kwargs:
            dicttrack.update({theparameter: kwargs[theparameter]})
    # default to parallel
    if "use_mp" not in dicttrack:
        dicttrack.update({"use_mp": True})

    # use radiation parameters to get the rf bucket
    pars = ring.radiation_parameters()

    # get energy bucket (deltabucket) max height.
    # S.Y.Lee, 4th Edition, Eqs 3.37, 3.38, Sect. II.2 Bucket area
    harmonnumber = ring.harmonic_number
    etac = pars.etac
    theenergy = pars.E0
    phis = pars.phi_s
    betar = ring.beta
    thevoltage = ring.get_rf_voltage()
    yfactor = numpy.sqrt(
        numpy.abs(numpy.cos(phis) - 0.5 * (numpy.pi - 2 * phis) * numpy.sin(phis))
    )
    deltabucket = (
        numpy.sqrt(
            2
            * thevoltage
            / (numpy.pi * betar**2 * theenergy * harmonnumber * numpy.abs(etac))
        )
        * yfactor
    )
    verboseprint(f"Bucket height {deltabucket}")

    # first guess
    es_ini = 0
    et_ini = deltabucket / 2
    if "eu_ini" in kwargs:
        eu_ini = kwargs["eu_ini"]
        verboseprint(f"Using the users max boundary {eu_ini}")
    else:
        eu_ini = deltabucket
        verboseprint("Using the bucket height as maximum boundary")

    if "orbit" in kwargs:
        orbit_s = kwargs["orbit"]
        verboseprint("Using the users orbit")
    else:
        _, orbit_s = ring.find_orbit(rps)
        verboseprint("Using the closed orbit")

    # start scan
    # deltaeu, unstable energy
    # deltaes, stable energy
    # deltaet, test energy
    etnp = numpy.ones((2, nrps))
    thesign = [-1, 1]  # first negative side, then positive side
    verbosesign = ["negative", "positive"]
    for i in [0, 1]:
        verboseprint(f"Search on the {verbosesign[i]} side")
        deltaet = thesign[i] * et_ini * numpy.ones(nrps)
        deltaeu = thesign[i] * eu_ini * numpy.ones(nrps)
        deltaes = thesign[i] * es_ini * numpy.ones(nrps)
        deltae = 1
        iteration = 0
        while iteration < 100 and (deltae > dptol):  # safety limit on iterations
            t00 = time.time()
            # plost is a mask, True for lost particles
            plost = multirefpts_track_islost(
                ring, rps, deltaet, orbit_s, nturns, add_offset, **dicttrack
            )
            deltaes[~plost] = deltaet[~plost]
            deltaeu[plost] = deltaet[plost]
            deltaet = (deltaes + deltaeu) / 2
            deltae = max(abs(deltaes - deltaeu))
            iteration = iteration + 1
            outmsg = (
                f"Iteration {iteration} in {verbosesign[i]} side",
                f" took {format(time.time()-t00):.3} s.",
                f" {deltae=}, {dptol=}",
            )
            verboseprint("".join(outmsg))
        etnp[i, :] = deltaet
    return etnp.T


def projectrefpts(ring, startrefpts, particles, **kwargs):
    """
    :py:fun:`projectrefpts` tracks particles from multiple reference
    points to a single end point.

    Usage:
      >>> projectrefpts(ring, startrefpts, particles)

    Parameters:
      ring: list of elements
      startrefpts: reference point to start tracking from.
      particles: (6,N,R,1) array, where N is the number of
        particles per reference point, and R is the number of
        reference points.

    Keyword arguments:
      endrefpt: end reference point. Default: end of last ring element
      use_mp: Default True. See :py:fun:`lattice_track`
      group: Default False. The starting point info is removed.
        All tracked particles are grouped together.
      verbose: prints additional info

    Returns:
      zout:     Tracked particle coordinates at the end ref. point.
                Default (6,N,R,1) array, where N is the number of
                  particles per R references points
                If the flag 'group' is used the output becomes a
                  (6,N*R,1,1) array with all particles.
      lostpart: Bool array, True when the particle is lost.
                Default (N,R).
                If the flag 'group' is used the output becomes a
                  (N*R) array
    """
    lenring = len(ring)
    rps = startrefpts
    nrps = len(rps)
    nparticles = 1
    if len(particles.shape) >= 2:
        nparticles = particles.shape[1]

    # verboseprint to check flag only once
    verbose = kwargs.pop("verbose", False)
    verboseprint = print if verbose else lambda *a, **k: None

    if "endrefpt" in kwargs:
        erps = kwargs["endrefpt"]
        verboseprint(f"Project particles to start of element index {erps}")
    else:
        erps = lenring
        verboseprint("Project to end point")

    groupparts = kwargs.pop("group", False)

    verboseprint(f"{nparticles=} per reference point")
    verboseprint(f"Number of reference points {nrps}")

    # default to parallel
    use_mp = kwargs.pop("use_mp", True)
    if nparticles == 1:
        use_mp = False

    zin = numpy.zeros((6, nparticles, nrps, 1))
    if groupparts:
        zout = numpy.zeros((6, nparticles * nrps, 1, 1))
        lostpart = numpy.ones((nparticles * nrps), dtype=bool)
    else:
        zout = numpy.zeros((6, nparticles, nrps, 1))
        lostpart = numpy.ones((nparticles, nrps), dtype=bool)

    # first, track the remaining portion of the ring
    zin = particles.copy()
    for i in range(nrps):
        ring_downstream = ring.rotate(rps[i])
        zaux = numpy.squeeze(zin[:, :, i, 0])
        verboseprint(
            f"Tracking {nparticles} particles on reference point {i+1} of {nrps}"
        )
        zoaux, _, dout1 = ring_downstream.track(
            zaux, nturns=1, refpts=erps - rps[i], losses=True, use_mp=use_mp
        )
        if groupparts:
            zout[:, nparticles * i : nparticles * (i + 1), 0, 0] = numpy.squeeze(zoaux)
            lostpart[nparticles * i : nparticles * (i + 1)] = dout1["loss_map"][
                "islost"
            ]
        else:
            zout[:, :, i, 0] = numpy.reshape(zoaux, (6, nparticles))
            lostpart[:, i] = dout1["loss_map"]["islost"]
    return zout, lostpart


def multirefpts_track_islost(
    ring: list,
    refpts: numpy.ndarray,
    energysetpt: float,
    orbit: numpy.ndarray,
    nturns: int,
    initcoord: numpy.ndarray,
    **dicttrack: dict[str, any],
) -> numpy.ndarray:
    """
    Returns a boolean array: tells whether the particle launched is lost.

    True means lost.
    """
    lenring = len(ring)
    rps = refpts
    nrps = len(rps)
    lostpart = numpy.ones((nrps), dtype=bool)
    zin = numpy.zeros((6, nrps))
    zout = numpy.zeros((6, nrps))
    issmall = 1e-6
    eps = numpy.finfo(float).eps
    istiny = 100 * eps

    # first, track the remaining portion of the ring
    for i in range(nrps):
        ring_downstream = ring.rotate(rps[i])  # [rp[i]:]
        zin[:, i] = orbit[i, :].copy()
        zin[0, :] = zin[0, :] + initcoord[0, :]
        zin[2, :] = zin[2, :] + initcoord[1, :]
        zin[4, :] = zin[4, :] + energysetpt
        zoaux, _, dout1 = ring_downstream.track(
            zin[:, i], nturns=1, refpts=lenring - rps[i], losses=True
        )
        lostpart[i] = dout1["loss_map"]["islost"][0]
        zout[:, i] = zoaux[:, 0, 0, 0]

    cntalive = nrps - sum(lostpart)
    aliveatringend = ~lostpart
    zinaliveaux = numpy.squeeze(zout[:, aliveatringend])
    zinalive = numpy.asfortranarray(zinaliveaux.copy())

    # second, track particles that have survived the ring to the end
    if cntalive == 1:
        _, _, dout2 = ring.track(zinalive, nturns=nturns, refpts=len(ring), losses=True)
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
        _, _, dout3 = ring.track(
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

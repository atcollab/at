"""momentum aperture alternative."""

from __future__ import annotations

import time
from typing import Dict

import numpy

from ..lattice.lattice_object import Lattice

__all__ = ["momaperture_project2start", "projectrefpts"]

# 2024may29 oblanco at ALBA CELLS. First working version
#           Based on the MATLAB implementation by Z.Marti at ALBA
#           See https://github.com/atcollab/at/pull/773


def momaperture_project2start(ring: Lattice, **kwargs: Dict[str, any]) -> numpy.ndarray:
    """
    :py:func:`momap_project2start` calculates the local momemtum aperture.

    It is a binary search of the negative and positive momentum thresholds
    of stability around the closed orbit.

    For a given momentum offset the particles are first tracked from every
    reference point to the end of the ring, and then, all or a set of particles
    with different 6D coordinates are tracked together.  The surviving particles
    continue the boundary search with a new energy step until the boundary is
    found when the step limit or the momentum limit is reached.

    Usage:
      >>> momaperture_project2start(ring)

    Parameters:
      ring: list of elements

    Keyword Arguments:
      refpts: Selects the location of coordinates output.
        See ":ref:`Selecting elements in a lattice <refpts>`"
      nturns: number of turns to be tracked. Default 1000
      dptol: momentum offset resolution. Default 1e-4
      dpuguess: maximum momentum boundary. Default: the rf bucket height.
      troffset: (2, N) offsets to be added to the transverse coordinates
        on the N reference points. Default 1e-5 m
      orbit: (N,6) offsets to be added on the N reference points.
        Default, the closed orbit
      verbose: print in the standard output additional info. Default False
      epsilon6d: float.
        If not passed, all particles are tracked.
        If epsilon6d is given, we track for nturns only particles
        having 6D coordinates different by epsilon6d

    Returns:
      dnp: (N,2) array with negative and positive stable momentum boundaries
        for the N reference points

    ..note::
      * This function could track in parallel. Set use_mp=True.
        Other arguments could be passed, check :py:func:`lattice_track`.
      * This function does a quick search, but, it is succeptible to miss
        islands of instability due to the varying energy step.
    """
    # verboseprint to check flag only once
    verbose = kwargs.pop("verbose", False)
    verboseprint = print if verbose else lambda *a, **k: None

    rps = kwargs.pop("refpts", ring.uint32_refpts(range(len(ring))))
    rps = ring.get_uint32_index(rps)
    nrps = len(rps)
    verboseprint(f"Using {nrps} reference points")

    nturns = kwargs.pop("nturns", 1000)
    verboseprint(f"Track over {nturns} turns")

    dptol = kwargs.pop("dptol", 1e-4)
    verboseprint(f"Momentum resolution {dptol}")

    # set transverse offsets
    if "troffset" in kwargs:
        add_offset = kwargs.pop("troffset")
        verboseprint("Add user offsets")
    else:
        dxy = 1e-5
        add_offset = numpy.tile(dxy, [2, nrps])
        verboseprint(f"Adding default transverse offsets {dxy} per plane")

    # set the minimum distance btw particles that makes them similar
    epsilon6d = kwargs.pop("epsilon6d", 0)

    # first guess
    if "dpuguess" in kwargs:
        eu_ini = kwargs.pop("dpuguess")
        verboseprint(f"Using the users max boundary {eu_ini}")
    else:
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
        eu_ini = deltabucket
        verboseprint("Using the bucket height as maximum boundary")
    es_ini = 0
    et_ini = eu_ini / 2

    if "orbit" in kwargs:
        orbit_s = kwargs.pop("orbit")
        verboseprint("Using the users orbit")
    else:
        _, orbit_s = ring.find_orbit(rps)
        verboseprint("Using the closed orbit")

    # start scan of:
    # unstable energy
    # stable energy
    # test energy
    etpos = et_ini * numpy.ones(nrps)
    eupos = eu_ini * numpy.ones(nrps)
    espos = es_ini * numpy.ones(nrps)
    etneg = -1 * et_ini * numpy.ones(nrps)
    euneg = -1 * eu_ini * numpy.ones(nrps)
    esneg = -1 * es_ini * numpy.ones(nrps)
    deltae = 1
    iteration = 0
    while iteration < 100 and (deltae > dptol):  # safety limit on iterations
        iteration = iteration + 1
        t00 = time.time()
        # plost is a mask, True for lost particles
        plost = multirefpts_track_islost(
            ring,
            rps,
            etpos,
            etneg,
            orbit_s,
            add_offset,
            nturns,
            epsilon6d,
            verbose,
            **kwargs,
        )
        # split in positive and negative side of energy offsets
        plostpos = plost[0::2]
        plostneg = plost[1::2]
        # split in stable (es), unstable (eu) and test (et) energy
        espos[~plostpos] = etpos[~plostpos]
        eupos[plostpos] = etpos[plostpos]
        etpos = (espos + eupos) / 2
        esneg[~plostneg] = etneg[~plostneg]
        euneg[plostneg] = etneg[plostneg]
        etneg = (esneg + euneg) / 2
        # define new energy step
        depos = max(abs(espos - eupos))
        deneg = max(abs(esneg - euneg))
        deltae = max(depos, deneg)
        outmsg = (
            f"Iteration {iteration}",
            f" took {format(time.time()-t00):.3} s.",
            f" dp_step={deltae}, dptol={dptol}",
        )
        verboseprint("".join(outmsg))
    return numpy.vstack([etneg, etpos]).T


def projectrefpts(
    ring: Lattice,
    startrefpts: numpy.ndarray,
    particles: numpy.ndarray,
    **kwargs: Dict[str, any],
) -> tuple:
    """
    :py:fun:`projectrefpts` tracks from multiple reference points.

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

    verboseprint(f"nparticles={nparticles} per reference point")
    verboseprint(f"Number of reference points {nrps}")

    # default to parallel
    use_mp = kwargs.pop("use_mp", True)
    if nparticles == 1:
        use_mp = False

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
            zout[:, nparticles * i : nparticles * (i + 1), 0, 0] = numpy.reshape(
                zoaux, (6, nparticles)
            )
            lostpart[nparticles * i : nparticles * (i + 1)] = dout1["loss_map"][
                "islost"
            ]
        else:
            zout[:, :, i, 0] = numpy.reshape(zoaux, (6, nparticles))
            lostpart[:, i] = dout1["loss_map"]["islost"]
    return zout, lostpart


def multirefpts_track_islost(
    ring: Lattice,
    refpts: numpy.ndarray,
    esetptpos: numpy.ndarray,
    esetptneg: numpy.ndarray,
    orbit: numpy.ndarray,
    initcoord: numpy.ndarray,
    nturns: float,
    epsilon6d: float,
    verbose: bool,
    **kwargs: Dict[str, any],
) -> numpy.ndarray:
    """
    Tell whether the particle launched is lost.

    Usage:
      >>> multirefpts_track_islost(ring, refpts, energysetpt, orbit, initcoord)

    Parameters:
      ring: list of elements

    Keyword Arguments:
      refpts: Selects the locations.
      esetptpos: positive energy set point for tracking.
      esetptneg: negative energy set point for tracking.
      orbit: (6,N) orbit to be added to the N refpts.
      initcoords: (2,N) hor. and ver. transverse offsets in m.
      nturns: number of turns to track
      epsilon6d: maximum value to consider particles as similar
      verbose: print info

    Returns:
      Lostpart: (2*N) bool array, for N reference points.
        True if the particle is lost.
    """
    verboseprint = print if verbose else lambda *a, **k: None

    # track positive and negative side at the same time
    nparticles = 2

    rps = refpts
    nrps = len(rps)
    lostpart = numpy.ones((nparticles * nrps), dtype=bool)
    zin = numpy.zeros((6, nparticles * nrps))
    zout = numpy.zeros((6, nparticles * nrps))
    tinyoffset = epsilon6d

    # project to the end of the ring
    erps = len(ring)

    # first, track the remaining portion of the ring
    zin[:, 0::2] = orbit.T.copy()
    zin[:, 1::2] = orbit.T.copy()
    zin[0, 0::2] = zin[0, 0::2] + initcoord[0, :]
    zin[0, 1::2] = zin[0, 1::2] + initcoord[0, :]
    zin[2, 0::2] = zin[2, 0::2] + initcoord[1, :]
    zin[2, 1::2] = zin[2, 1::2] + initcoord[1, :]
    zin[4, 0::2] = zin[4, 0::2] + esetptpos
    zin[4, 1::2] = zin[4, 1::2] + esetptneg
    for i in range(nrps):
        ring_downstream = ring.rotate(rps[i])
        zaux = zin[:, (nparticles * i) : (nparticles * (i + 1))]
        # verboseprint(
        #    f"Tracking {nparticles} particles on reference point {i+1} of {nrps}"
        # )
        zoaux, _, doutaux = ring_downstream.track(
            zaux, nturns=1, refpts=erps - rps[i], losses=True, **kwargs
        )
        zout[:, (nparticles * i) : (nparticles * (i + 1))] = numpy.reshape(
            zoaux, (6, nparticles)
        )
        lostpart[(nparticles * i) : (nparticles * (i + 1))] = doutaux["loss_map"][
            "islost"
        ]

    cntalive = len(lostpart) - sum(lostpart)
    zinaliveaux = zout[:, ~lostpart]
    zinalive_at_ring_end = numpy.asfortranarray(zinaliveaux.copy())

    # second, use the particles that have survived the ring to the end
    # filter them if necessary, and track them
    shapealiveatend = numpy.shape(zinalive_at_ring_end)
    trackonly_mask = numpy.ones(shapealiveatend[1], dtype=bool)
    similarparticles_index = numpy.array([])
    particles_were_filtered = False
    if epsilon6d != 0 and cntalive > 1:
        particles_were_filtered = True
        # search for non numerically similar particles
        closenessmatrix = numpy.zeros((cntalive, cntalive), dtype=bool)
        for i in range(cntalive):
            closenessmatrix[i, :] = [
                numpy.allclose(
                    zinalive_at_ring_end[:, j],
                    zinalive_at_ring_end[:, i],
                    atol=tinyoffset,
                )
                for j in range(cntalive)
            ]
        _, rowidx = numpy.indices((cntalive, cntalive))
        maxidx = numpy.max(rowidx * closenessmatrix, 1)
        trackonly_mask, _, similarparticles_index = numpy.unique(
            maxidx, return_index=True, return_inverse=True
        )
        outmsg = (
            "Speed up when discarding similar particles, ",
            f"{100*len(trackonly_mask)/cntalive:.3f}%",
        )
        verboseprint("".join(outmsg))
    # track
    # dummy track to later reuse the ring
    zaux2 = numpy.array([initcoord[0, 0], 0, initcoord[1, 0], 0, 1e-6, 0]).T
    ring.track(zaux2, nturns=1)
    # track non-numerically similar particles
    _, _, dout_multiturn = ring.track(
        zinalive_at_ring_end[:, trackonly_mask],
        nturns=nturns,
        keep_lattice=True,
        losses=True,
        **kwargs,
    )
    if particles_were_filtered:
        lostpaux = dout_multiturn["loss_map"]["islost"][similarparticles_index]
    else:
        lostpaux = dout_multiturn["loss_map"]["islost"]
    lostpart[~lostpart] = lostpaux

    return lostpart

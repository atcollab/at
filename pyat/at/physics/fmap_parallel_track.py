"""
Frequency analysis (FMAP) using PyNAFF lib and pyat parallel tracking (patpass)
"""

# orblancog
# generates the frequency and diffusion map for a given ring
# 2023jan16 tracking is parallel (patpass), NAFF analysis is serial
# 2022jun07 serial version

import at
from   at.tracking import patpass
import numpy                      as np
import PyNAFF

__all__ = ['fmap_parallel_track']

def fmap_parallel_track(ring, \
        coords = [-10,10,-3,3], \
        step = [0.05,0.05], \
        turns = 512, \
        ncpu = 30, \
        co = False, \
        add_offset6D = np.zeros((6,1)), \
        verbose = False, \
        ):
    """
    This function calculates the norm of the transverse tune variation per turn
    for a particle tracked along a ring with a set of offsets in the initial
    coordinates.

    It returns a numpy array containing 5 columns,
        [xoffset, yoffset, nux, nuy, log10( sqrt(dnux*dnux + dnuy*dnuy)/tns )]
    for every tracked particle that survives 2*tns turns.
    Particles lost before 2*tns turns are ignored.

    The transverse offsets are given inside a rectangular coordinate window
        coords=[xmin,xmax,ymin,ymax]
    in milimeters, where the window is divided in steps of step=[xstep,ystep].
    For each yoffset, particle tracking over the whole set of xoffsets is done
    in parallel, therefore, a number of cpus (ncpu) above the number of xsteps
    will not reduce the calculation time.

    The frequency analysis uses the library PyNAFF (not parallelized).

    The closed orbit (co) could be calculated and added to the
    initial particle offset of every particle by setting co=True.
    Additionally, a six by one numpy array (add_offset6D) could be used to
    arbitrarily offset the initial coordinates of every particle.

    Parameters :
      ring:     a valid pyat ring
    Optional:
      coords:   default [-10,10,-3,3] in mm
      step:     default [0.05,0.05]   in mm
      turns:    default 512
      ncpu:     max. number of processors in parallel tracking patpass
      co:       default false
      add_offset6D: default np.zeros((6,1))
      verbose:  prints additional info
    Returns:
      tune_and_nudiff_array: numpy array with columns
                             [xcoor, ycoor, nux, ny, log10(sqrt(sum(dnu**2)/turns)) ]

    WARNING : points with NaN tracking results or non-defined frequency
              in x,y are ignored.
    WARNING : if co flag is used, the closed orbit is added using
                  at.physics.find_orbit
              however, it is not returned by this function.
    """

    if co:
        # closed orbit values are not returned. It seems not necessary here
        addco = at.physics.find_orbit(ring);
        print(f'Closed orbit:\t{addco}')
    else:
        addco = [np.zeros((6,1))]

    # simple verbose to check flag only once
    if verbose:
        def verboseprint(*args):
            # print each argument separately
            for arg in args: print arg, ; print
    else:
        verboseprint = lambda *a: None

    # tns is the variable used in the frequency analysis
    # turns is the input variable from user
    # nturns is twice tns in order to get the tune in the first and second
    #     part of the tracking
    tns = turns;
    nturns = 2*tns;

    # scale the input coordinates
    xscale = 1e-3;
    yscale = 1e-3;

    # returned array
    tune_and_nudiff_array = np.empty([])

    # define rectangle and x,y step size
    xmin  = coords[0]
    xmax  = coords[1]
    ymin  = coords[2]
    ymax  = coords[3]
    xstep = step[0]
    ystep = step[1]

    # get the intervals
    ixarray    = np.arange(xmin, xmax+1e-6, xstep)
    lenixarray = len(ixarray);
    iyarray    = np.arange(ymin, ymax+1e-6, ystep)
    leniyarray = len(iyarray);

    print("Start tracking and NAFF analysis")

    # Forcing a maximum value of cpus
    # I have no idea how to confirm the actual used value from patpass
    # at.DConstant.patpass_poolsize = ncpu;
    print(f' POOL size : {ncpu}')

    # tracking in parallel multiple x coordinates with the same y coordinate
    for iy,iy_index in zip(iyarray, range(leniyarray)):
          print(f'Tracked particles {abs(-100.0*iy_index/leniyarray):.1f} %, with a max. cpu POOL size of {ncpu}');
          print("y =",iy)
          #z01 = np.array([ix*xscale+1e-9, 0,  0, 0,  0, 0])
          # add 1 nm to tracking to avoid zeros in array for the ideal lattice
          z0 = np.zeros((6,lenixarray))
          #z0 = z0 + add_offset6D + addco[0,:];
          #z0[0,:] = z0[0,:] + xscale*ixarray + 1e-9;
          #z0[2,:] = z0[2,:] + yscale*iy      + 1e-9;
          z0[0,:] = xscale*ixarray + add_offset6D[0] + addco[0][0] + 1e-9 ;
          z0[1,:] =                + add_offset6D[1] + addco[0][1];
          z0[2,:] = yscale*iy      + add_offset6D[2] + addco[0][2] + 1e-9 ;
          z0[3,:] =                + add_offset6D[3] + addco[0][3];
          z0[4,:] =                + add_offset6D[4] + addco[0][4];
          z0[5,:] =                + add_offset6D[5] + addco[0][5];
          zOUT = patpass(ring, z0, nturns, pool_size=ncpu);

          # start of serial frequency analysis
          for ix_index in range(lenixarray): # cycle over the track results
              # check if nan in arrays
              array_sum = np.sum(zOUT[:,ix_index,0]);
              array_has_nan = np.isnan(array_sum)
              if array_has_nan:
                  verboseprint("array has nan")
                  continue


              # get one valid particle
              z1 = zOUT[:,ix_index,0];

              # remove mean values
              # get the first turn in x
              xfirst     = z1[0, 0:tns];
              xfirst     = xfirst - np.mean(xfirst);
              pxfirst    = z1[1, 0:tns];
              pxfirst    = pxfirst - np.mean(pxfirst);
              xfirstpart = xfirst + 1j*pxfirst;
              # get the last turns in x
              xlast      = z1[0, tns:2*tns];
              xlast      = xlast - np.mean(xlast);
              pxlast     = z1[1, tns:2*tns];
              pxlast     = pxlast - np.mean(pxlast);
              xlastpart  = xlast + 1j*pxlast;

              # get the first turn in y
              yfirst     = z1[2, 0:tns];
              yfirst     = yfirst - np.mean(yfirst);
              pyfirst    = z1[3, 0:tns];
              pyfirst    = pyfirst - np.mean(pyfirst);
              yfirstpart = yfirst + 1j*pyfirst;
              # get the last turns in y
              ylast      = z1[2, tns:2*tns];
              ylast      = ylast - np.mean(ylast);
              pylast     = z1[3, tns:2*tns];
              pylast     = pylast - np.mean(pylast);
              ylastpart  = ylast + 1j*pylast;

              # calc frequency from array,
              # jump the cycle is no frequency is found
              xfreqfirst = PyNAFF.naff(xfirstpart,tns,1,0,False)
              if len(xfreqfirst) == 0:
                  verboseprint("  No frequency");
                  continue;
              xfreqlast  = PyNAFF.naff(xlastpart,tns,1,0,False)
              if len(xfreqlast) == 0:
                  verboseprint("  No frequency");
                  continue;
              yfreqfirst = PyNAFF.naff(yfirstpart,tns,1,0,False)
              if len(yfreqfirst) == 0:
                  verboseprint("  No frequency");
                  continue;
              yfreqlast  = PyNAFF.naff(ylastpart,tns,1,0,False)
              if len(yfreqlast) == 0:
                  verboseprint("  No frequency");
                  continue;
              verboseprint("NAFF results", \
                           "\nH freq. first part =\t" ,xfreqfirst[0][1], \
                           "\nH freq. last part =\t", xfreqlast[0][1], \
                           "\nV freq. first part =\t", yfreqfirst[0][1], \
                           "\nV freq. last part =\t", yfreqlast[0][1])

              # metric
              xdiff = xfreqlast[0][1] - xfreqfirst[0][1];
              ydiff = yfreqlast[0][1] - yfreqfirst[0][1];
              nudiff = 0.5*np.log10((xdiff*xdiff + ydiff*ydiff)/tns)
              # min max diff
              if  nudiff >  -2:
                  nudiff =  -2;
              if  nudiff < -10:
                  nudiff = -10;
              # save diff
              tune_and_nudiff_array = np.append( \
                      tune_and_nudiff_array, \
                      ixarray[ix_index],iy, \
                      [xfreqfirst[0][1], yfreqfirst[0][1], \
                      nudiff]);

    # first element is garbage
    tune_and_nudiff_array = np.delete(tune_and_nudiff_array,0);
    ## reshape for plots and output files
    tune_and_nudiff_array = tune_and_nudiff_array.reshape(-1,5);

    return tune_and_nudiff_array
# the end


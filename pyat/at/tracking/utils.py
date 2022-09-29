"""
Utility functions for tracking simulations
"""
import numpy


def get_bunches(r_in, nbunch, selected_bunches=None):
    """Function to get the bunches particles 6D coordinates
    
    Parameters:
      r_in: 6 x n_particles Fortran-ordered numpy array.
      nbunch: integer, total number of bunches
      selected_bunches: integer or array of integers, index
      of the selected bunches
      
    Returns:
      List of ndarray containing the 6 x n particles coordinates
      of the selected bunches     
    """
    if selected_bunches is None:
        selected_bunches = numpy.arange(nbunch)
    else:
        selected_bunches = numpy.atleast_1d(selected_bunches)
    bunches = [r_in.T[ib::nbunch].T for ib in selected_bunches]
    return bunches
    
    
def get_bunches_std_mean(r_in, nbunch, selected_bunches=None):
    """Function to get the bunches standard deviation and center
    of mass
            
    Parameters:
      r_in: 6 x n_particles Fortran-ordered numpy array.
      nbunch: integer, total number of bunches
      selected_bunches: integer or array of integers, index
      of the selected bunches
      
    Returns:
      Lists of ndarray containing the 6D standard deviation
      and center of mass (std, mean) 
    """
    bunches = get_bunches(r_in, nbunch, selected_bunches)
    std = [numpy.nanstd(b, axis=1) for b in bunches]
    mean = [numpy.nanmean(b, axis=1) for b in bunches]
    return std, mean 
    
        



from warnings import warn
from .atpass import atpass as _atpass, elempass as _elempass


__all__ = ['atpass', 'elempass']


# noinspection PyIncorrectDocstring
def atpass(*args, **kwargs):
    """
    atpass(line, r_in, nturns, refpts=[], reuse=False, omp_num_threads=0)

    Track input particles *r_in* along line for *nturns* turns.

    Parameters:
        line (Sequence[Element]): list of elements
        r_in:                   6 x n_particles Fortran-ordered numpy array.
          On return, rin contains the final coordinates of the particles
        nturns (int):           number of turns to be tracked

    Keyword arguments:
        refpts (Uint32_refs):   numpy array of indices of elements where
          output is desired:

          * 0 means entrance of the first element
          * len(line) means end of the last element

        energy:                 nominal energy [eV]
        rest_energy:            rest_energy of the particle [eV]
        charge:                 particle charge [elementary charge]
        reuse (bool):           if True, use previously cached description
          of the lattice.
        omp_num_threads (int):  number of OpenMP threads
          (default 0: automatic)
        losses (bool):          if True, process losses

    Returns:
        r_out:  6 x n_particles x n_refpts x n_turns Fortran-ordered
          numpy array of particle coordinates

    :meta private:
    """
    warn(UserWarning("The public interface for tracking is 'lattice_pass'"))
    return _atpass(*args, **kwargs)


# noinspection PyIncorrectDocstring
def elempass(*args, **kwargs):
    """elempass(element, r_in)

    Track input particles *r_in* through a single element.

    Parameters:
        element (Element):  AT element
        rin:                6 x n_particles Fortran-ordered numpy array.
          On return, rin contains the final coordinates of the particles

    Keyword arguments:
        energy:             nominal energy [eV]
        rest_energy:        rest_energy of the particle [eV]
        charge:             particle charge [elementary charge]

    :meta private:
    """
    warn(UserWarning("The public interface for tracking is 'element_pass'"))
    return _elempass(*args, **kwargs)

from warnings import warn
from math import pi
import numpy
from scipy.constants import physical_constants as cst
import at       # for AtWarning, AtError
from ..lattice import elements

__all__ = ['get_energy', 'get_voltage', 'get_energy_loss']

TWO_PI_ERROR = 1.E-4


def get_energy(ring):
    """
    energy, periodicity = get_energy(ring) Get the ring energy
    get_energy looks for the machine energy in:
        1) the 1st 'RingParam' element
        2) the 1st 'RFCavity' element

    PARAMETERS
        ring  Ring structure

    OUPUT
        energy          Ring energy in eV
        periodicity     Number of periods to make 2pi bending
    """
    params = [elem for elem in ring if isinstance(elem, elements.RingParam)]

    if len(params) > 0:
        energy = params[0].Energy
        periodicity = params[0].Periodicity
    else:
        cavities = [elem for elem in ring if isinstance(elem, elements.RFCavity)]
        if len(cavities) > 0:
            energy = cavities[0].Energy
            theta = [elem.BendingAngle for elem in ring if isinstance(elem, elements.Dipole)]
            try:
                nbp = 2.0 * pi / sum(theta)
            except ZeroDivisionError:
                warn(at.AtWarning('No bending in the cell, set "Periodicity" to 1'))
                periodicity = 1
            else:
                periodicity = int(round(nbp))
                if abs(periodicity - nbp) > TWO_PI_ERROR:
                    warn(at.AtWarning('non-integer number of cells: {0} -> {1}'.format(nbp, periodicity)))
        else:
            raise at.AtError('Energy not defined (searched in "RingParam" and "Cavities")')

    return energy, periodicity


def get_voltage(ring):
    """
    voltage, harm_number = get_voltage(ring) Get the accelerating voltage

    PARAMETERS
        ring  Ring structure

    OUPUT
        voltage    Total RF voltage
        harmnumber Harmonic number
    """
    cavities = [elem for elem in ring if isinstance(elem, elements.RFCavity)]
    energy, periodicity = get_energy(ring)

    if len(cavities) > 0:
        voltage = periodicity * sum([elem.Voltage for elem in cavities])
        harm_number = periodicity * cavities[0].HarmNumber
    else:
        voltage = float('nan')
        harm_number = None

    return voltage, harm_number


def get_energy_loss(ring):
    """
    energy_loss = energy_loss(ring) Get the energy loss per turn

    PARAMETERS
        ring  Ring structure

    OUPUT
        energy_loss     energy loss per turn in eV

    Losses = Cgamma / 2pi * EGeV^4 * I2
    """
    lenthe = numpy.array([(elem.Length, elem.BendingAngle) for elem in ring if isinstance(elem, elements.Dipole)])
    lendp = lenthe[:, 0]
    theta = lenthe[:, 1]
    energy, periodicity = get_energy(ring)

    e_radius = cst['classical electron radius'][0]
    e_mass = cst['electron mass energy equivalent in MeV'][0]
    cgamma = 4.0E9*pi*e_radius/3.0/pow(e_mass, 3)

    I2 = periodicity*(numpy.sum(theta*theta/lendp))
    e_loss = cgamma/2.0/pi*pow(energy*1.0E-9, 4) * I2 * 1.e9
    return e_loss

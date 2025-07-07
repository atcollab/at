"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""

from __future__ import annotations

__all__ = [
    "LongtMotion",
    "_DictLongtMotion",
    "_Radiative",
    "Radiative",
    "Collective",
]

import abc
from abc import ABC
from copy import deepcopy


class LongtMotion(ABC):
    """Abstract Base class for all Element classes whose instances may modify
    the particle momentum

    Allows identifying elements potentially inducing longitudinal motion.

    Subclasses of :py:class:`LongtMotion` must provide two methods for
    enabling longitudinal motion:

    * ``_get_longt_motion(self)`` must return the activation state,
    * ``set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs)``
      must enable or disable longitudinal motion.
    """

    @abc.abstractmethod
    def _get_longt_motion(self):
        return False

    # noinspection PyShadowingNames
    @abc.abstractmethod
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        """Enable/Disable longitudinal motion

        Parameters:
            enable:     :py:obj:`True`: for enabling, :py:obj:`False` for
              disabling
            new_pass:   New PassMethod:

              * :py:obj:`None`: makes no change,
              * ``'auto'``: Uses the default conversion,
              * Anything else is used as the new PassMethod.
            copy:       If True, returns a modified copy of the element,
              otherwise modifies the element in-place
        """
        # noinspection PyUnresolvedReferences
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if copy:
            newelem = deepcopy(self)
            newelem.PassMethod = new_pass
            return newelem
        # noinspection PyAttributeOutsideInit
        self.PassMethod = new_pass
        return None


# noinspection PyUnresolvedReferences
class _DictLongtMotion(LongtMotion):
    # noinspection PyShadowingNames
    """Mixin class for elements implementing a 'default_pass' class attribute

    :py:class:`DictLongtMotion` provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod according
      to the ``default_pass`` dictionary.
    * a :py:obj:`.longt_motion` property set to :py:obj:`True` when the
      PassMethod is ``default_pass[True]``

    The class must have a ``default_pass`` class attribute, a dictionary
    such that:

    * ``default_pass[False]`` is the PassMethod when radiation is turned
      OFF,
    * ``default_pass[True]`` is the default PassMethod when radiation is
      turned ON.

    The :py:class:`DictLongtMotion` class must be set as the first base class.

    Example:

        >>> class QuantumDiffusion(_DictLongtMotion, Element):
        ...     default_pass = {False: "IdentityPass", True: "QuantDiffPass"}

        Defines a class such that :py:meth:`set_longt_motion` will select
        ``'IdentityPass'`` or ``'IdentityPass'``.
    """

    def _get_longt_motion(self):
        return self.PassMethod != self.default_pass[False]

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == "auto":
            new_pass = self.default_pass[enable]
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)


# noinspection PyUnresolvedReferences
class _Radiative(LongtMotion):
    # noinspection PyShadowingNames
    r"""Mixin class for radiating elements

    :py:class:`_Radiative` implements the mechanism for converting the pass
    methods of radiating elements. It provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod
      according to the following rule:

      * ``enable == True``: replace "\*Pass" by "\*RadPass"
      * ``enable == False``: replace "\*RadPass" by "\*Pass"
    * a :py:obj:`.longt_motion` property set to true when the PassMethod
      ends with "RadPass"

    The :py:class:`_Radiative` class must be set as the first base class.

    Example:
        >>> class Multipole(_Radiative, LongElement, ThinMultipole):

        Defines a class where :py:meth:`set_longt_motion` will convert the
        PassMethod according to the \*Pass or \*RadPass suffix.
    """

    def _get_longt_motion(self):
        return self.PassMethod.endswith(("RadPass", "QuantPass"))

    def _autopass(self, enable):
        if enable:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return "".join((root[:-4], "RadPass"))
        elif self.longt_motion:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return root
        else:
            return None

    # noinspection PyTypeChecker,PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        if new_pass == "auto":
            new_pass = self._autopass(enable)
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if enable:

            def setpass(el):
                el.PassMethod = new_pass
                el.Energy = kwargs["energy"]

        else:

            def setpass(el):
                el.PassMethod = new_pass
                try:
                    del el.Energy
                except AttributeError:
                    pass

        if copy:
            newelem = deepcopy(self)
            setpass(newelem)
            return newelem
        setpass(self)
        return None


class Radiative(_Radiative):
    # noinspection PyUnresolvedReferences
    r"""Mixin class for default radiating elements (:py:class:`.Dipole`,
    :py:class:`.Quadrupole`, :py:class:`.Wiggler`)

    :py:class:`Radiative` is a base class for the subset of radiative elements
    considered as the ones to be turned on by default: :py:class:`.Dipole`,
    :py:class:`.Quadrupole` and :py:class:`.Wiggler`, excluding the higher
    order multipoles.

    :py:class:`Radiative` inherits from :py:class:`_Radiative` and does not
    add any new functionality. Its purpose is to identify the default set of
    radiating elements.

    Example:
        >>> class Dipole(Radiative, Multipole):

        Defines a class belonging to the default radiating elements. It
        converts the PassMethod according to the "\*Pass" or "\*RadPass"
        suffix.
    """

    pass


class Collective(_DictLongtMotion):
    # noinspection PyAbstractClass,PyUnresolvedReferences
    """Mixin class for elements representing collective effects

    Derived classes will automatically set the
    :py:attr:`~Element.is_collective` property when the element is active.

    The class must have a ``default_pass`` class attribute, a dictionary such
    that:

    * ``default_pass[False]`` is the PassMethod when collective effects
      are turned OFF,
    * ``default_pass[True]`` is the default PassMethod when collective effects
      are turned ON.

    The :py:class:`Collective` class must be set as the first base class.

    Example:
        >>> class WakeElement(Collective, Element):
        ...     default_pass = {False: "IdentityPass", True: "WakeFieldPass"}

        Defines a class where the :py:attr:`~Element.is_collective` property is
        handled
    """

    def _get_collective(self):
        # noinspection PyUnresolvedReferences
        return self.PassMethod != self.default_pass[False]

    @abc.abstractmethod
    def clear_history(self):
        pass

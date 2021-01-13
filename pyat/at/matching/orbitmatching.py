from at.physics import find_orbit4
from at.matching import ElementConstraints


class Orbit4Constraints(ElementConstraints):
    """Container for orbit constraints:
      at.find_orbit4 is called once before the evaluation of all constraints

      Example:
          cnstrs = Orbit4Constraints(ring, dp=0.01)

          # Add a bump (x=-0.004, x'=0) constraint at location ref_inj
          cnstrs.add([-0.004, 0.0], refpts=ref_inj, order=slice(2))
    """
    def __init__(self, ring, *args, **kwargs):
        """Build a Orbit4Constraints container

        KEYWORDS
        dp=0.0              momentum deviation.
        orbit=None          Initial trajectory for transfer line: (6,) array
        """
        super(Orbit4Constraints, self).__init__(ring, *args, **kwargs)

    def add(self, target, refpts=None, order=None, name=None, **kwargs):
        """Add a target to the Orbit4Constraints container

        PARAMETERS
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint. Several locations may be
                          given to apply the same constraint at several points.
            order=None    index in the orbit vector. If None, the full orbit
                          is used. Example:
                            order=0         # x
                            order=2         # z
                            order=slice(4)  # x, x', z, z'
            name='orbit4' name of the constraint.
            weight=1.0    weight factor: the residual is (value-target)/weight.
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """

        if name is None:                # Generate the constraint name
            name = 'orbit4' if order is None else 'orbit4_{}'.format(order)

        fun = self._arrayaccess(order)
        super(Orbit4Constraints, self).add(fun, target, refpts, name=name,
                                           **kwargs)

    def compute(self, ring, *args, **kwargs):
        """Orbit computation before evaluation of all constraints"""
        orbit0, orbit = find_orbit4(ring, refpts=self.refpts, **kwargs)
        return (orbit[ref[self.refpts]].T for ref in self.refs), ()

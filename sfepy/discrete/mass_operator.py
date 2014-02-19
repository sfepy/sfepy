from sfepy.base.base import Struct
from sfepy.solvers import Solver

class MassOperator(Struct):
    """
    Encapsulation of action and inverse action of a mass matrix operator
    :math:`M`.
    """

    def __init__(self, problem, options):
        self.mtx_mass = problem.evaluate(options.mass, mode='weak',
                                         auto_init=True, dw_mode='matrix')

        if options.lumped:
            raise NotImplementedError

        else:
            # Initialize solvers (and possibly presolve the matrix).
            self.ls = Solver.any_from_conf(problem.ls_conf, mtx=self.mtx_mass,
                                           presolve=True)

    def action(self, vec):
        """
        Action of mass matrix operator on a vector: :math:`M x`.
        """
        return self.mtx_mass * vec

    def inverse_action(self, vec):
        """
        Inverse action of mass matrix operator on a vector: :math:`M^{-1} x`.
        """
        return self.ls(vec)

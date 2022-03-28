#  hIPPYlib-MUQ interface for large-scale Bayesian inverse problems
#  Copyright (c) 2019-2020, The University of Texas at Austin,
#  University of California--Merced, Washington University in St. Louis,
#  The United States Army Corps of Engineers, Massachusetts Institute of Technology

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
import dolfin as dl
import hippylib as hp
from minimization import NewtonFwdSolver


class NonlinearPPossionForm(object):
    def __init__(self, p, f, ds_n):
        """
        :param p: the exponent
        :param f: the forcing term
        :param ds_n: measure of the boundary on which the Neumann condition applies
        """
        self.p = p
        self.f = f
        self.ds_n = ds_n

        # Regularization term
        self.eps = dl.Constant(1e-8)

    def energy_functional(self, u, m):
        """Energy functional

        :param u: state variable
        :param m: parameter variable (here the Neumann data)
        :returns: energy functional
        """
        grad_u = dl.nabla_grad(u)
        etah = dl.inner(grad_u, grad_u) + self.eps

        if self.f is None:
            return 1.0 / self.p * etah ** (0.5 * self.p) * dl.dx - u * m * self.ds_n
        else:
            return (
                1.0 / self.p * etah ** (0.5 * self.p) * dl.dx
                + self.f * u * dl.dx
                - u * m * self.ds_n
            )

    def variational_form(self, u, m, p):
        """Variational form

        :param u: state variable
        :param m: parameter variable
        :param p: adjoint variable
        :returns: variational form
        """
        return dl.derivative(self.energy_functional(u, m), u, p)


class EnergyFunctionalPDEVariationalProblem(hp.PDEVariationalProblem):
    def __init__(self, Vh, energyform, bc, bc0):
        hp.PDEVariationalProblem.__init__(
            self, Vh, energyform.variational_form, bc, bc0
        )

        self.energy_fun = energyform.energy_functional
        self.fwd_solver = NewtonFwdSolver()
        self.qoi = None
        self.cal_qoi = False

    def solveFwd(self, state, x):
        """Solve the nonlinear forward problem using Newton's method

        :param state: state variable
        :param x: x[0] = state, x[1] = parameter, x[2] = adjoint variables
        """
        if self.solver is None:
            self.solver = self._createLUSolver()

        u = hp.vector2Function(x[hp.STATE], self.Vh[hp.STATE])
        m = hp.vector2Function(x[hp.PARAMETER], self.Vh[hp.PARAMETER])

        F = self.energy_fun(u, m)

        uvec, reason = self.fwd_solver.solve(F, u, self.bc, self.bc0)

        if not self.fwd_solver.converged:
            print("Newton did not converged", reason)

        state.zero()
        state.axpy(1.0, uvec.vector())

        if self.cal_qoi:
            self.qoi.update_tracer(state)

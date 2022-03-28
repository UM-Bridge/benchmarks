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

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
The basic part of this code is taken from
https://github.com/hippylib/hippylib/blob/master/applications/poisson/model_subsurf.py.

Input values of this script should be defined in "ppoisson_box.yaml"
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(""))))

import math
import yaml
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import dolfin as dl
import hippylib as hp
import muq.Modeling as mm
import muq.SamplingAlgorithms as ms
import hippylib2muq as hm

from nonlinearPPoissonProblem import *


def true_model(prior):
    """
    Define true parameter field.

    In this example, we sample from the prior and take it as the true parameter
    field.
    """
    noise = dl.Vector()
    prior.init_vector(noise, "noise")
    hp.parRandom.normal(1.0, noise)
    mtrue = dl.Vector()
    prior.init_vector(mtrue, 0)
    prior.sample(noise, mtrue)
    return mtrue


def export2XDMF(x, Vh, fid):
    fid.parameters["functions_share_mesh"] = True
    fid.parameters["rewrite_function_mesh"] = False

    fun = hp.vector2Function(x, Vh)
    fid.write(fun, 0)


class BottomBoundary(dl.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and dl.near(x[2], 0)


class SideBoundary(dl.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (
            dl.near(x[0], 0)
            or dl.near(x[0], Length)
            or dl.near(x[1], 0)
            or dl.near(x[1], Width)
        )


class TopBoundary(dl.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and dl.near(x[2], Height)


def generate_starting():
    """
    Generate an initial parameter sample from the Laplace posterior for the MUQ
    MCMC simulation
    """
    noise = dl.Vector()
    nu.init_vector(noise, "noise")
    hp.parRandom.normal(1.0, noise)
    pr_s = model.generate_vector(hp.PARAMETER)
    post_s = model.generate_vector(hp.PARAMETER)
    nu.sample(noise, pr_s, post_s, add_mean=True)
    x0 = hm.dlVector2npArray(post_s)
    return x0


def data_file(action, target=None, data=None):
    """
    Read or write the observations.

    :param action: "w" is to write the date to the file named "data.h5" and "r"
                   is to read the data from "data.h5"
    :param target: the location of the observation data
    :param data: the observation data
    """
    f = h5py.File("data.h5", action)
    if action == "w":
        f["/target"] = target
        f["/data"] = data

        f.close()
        return

    elif action == "r":
        target = f["/target"][...]
        data = f["/data"][...]

        f.close()

        return target, data


class ExtractBottomData:
    def __init__(self, mesh, Vh):
        bmesh = dl.BoundaryMesh(mesh, "exterior")
        bmarker = dl.MeshFunction("size_t", mesh, bmesh.topology().dim())
        for c in dl.cells(bmesh):
            if math.isclose(c.midpoint().z(), 0):
                bmarker[c] = 1

        smesh = dl.SubMesh(bmesh, bmarker, 1)

        self.vertex_s2b = smesh.data().array("parent_vertex_indices", 0)
        self.vertex_b2p = bmesh.entity_map(0).array()
        self.vertex2dof = dl.vertex_to_dof_map(Vh)
        self.coordinates = smesh.coordinates()

        self.tria = tri.Triangulation(
            self.coordinates[:, 0], self.coordinates[:, 1], smesh.cells()
        )

    def get_dim(self):
        return self.coordinates.shape[0]

    def get_bottom_data(self, arr):
        return arr[self.vertex2dof[self.vertex_b2p[self.vertex_s2b]]]

    def plot_array(self, arr, vmin=None, vmax=None, cmap=None, fname=None):
        val = arr[self.vertex2dof[self.vertex_b2p[self.vertex_s2b]]]

        if vmax is None:
            vmax = np.max(val)
        if vmin is None:
            vmin = np.min(val)

        plt.tripcolor(self.tria, val, shading="gouraud", vmin=vmin, vmax=vmax)
        if cmap:
            plt.set_cmap(cmap)

        plt.axis("off")
        plt.gca().set_aspect("equal")

        if fname:
            plt.savefig(fname, dpi=100, bbox_inches="tight", pad_inches=0)

        plt.show()


class TracerSideFlux:
    def __init__(self, ds, p, n):
        self.n = dl.FacetNormal(mesh)
        self.ds = ds
        self.p = p

        self.tracer = hp.QoiTracer(n)
        self.ct = 0

    def form(self, u):
        grad_u = dl.nabla_grad(u)
        etah = dl.inner(grad_u, grad_u)

        return etah ** (0.5 * (self.p - 2)) * dl.dot(grad_u, self.n) * self.ds

    def eval(self, u):
        uf = hp.vector2Function(u, Vh[hp.STATE])
        return dl.assemble(self.form(uf))

    def update_tracer(self, state):
        y = self.eval(state)
        self.tracer.append(self.ct, y)
        self.ct += 1


def paramcoord2eigencoord(V, B, x):
    """
    Projection a parameter vector to eigenvector.

    y = V^T * B * x

    :param V multivector: eigenvectors
    :param operator: the right-hand side operator in the generalized eig problem
    :param x np.array: parameter data
    """
    # convert np.array to multivector
    nvec = 1
    Xvecs = hp.MultiVector(pde.generate_parameter(), nvec)
    hm.npArray2dlVector(x, Xvecs[0])

    # multipy B
    BX = hp.MultiVector(Xvecs[0], nvec)
    hp.MatMvMult(B, Xvecs, BX)
    VtBX = BX.dot_mv(V)

    return VtBX.transpose()


if __name__ == "__main__":
    with open("ppoisson_box.yaml") as fid:
        inargs = yaml.full_load(fid)

    sep = "\n" + "#" * 80 + "\n"

    #
    #  Set up the mesh and finite element function spaces
    #
    ndim = 3
    Length = 1.0
    Width = Length
    Height = 0.05

    nx = inargs["nelement"][0]
    ny = nx
    nz = inargs["nelement"][1]

    mesh = dl.BoxMesh(dl.Point(0, 0, 0), dl.Point(Length, Width, Height), nx, ny, nz)
    bottom = BottomBoundary()
    side = SideBoundary()
    top = TopBoundary()

    Vh1 = dl.FunctionSpace(mesh, "Lagrange", 1)
    Vh2 = dl.FunctionSpace(mesh, "Lagrange", 1)
    Vh = [Vh2, Vh1, Vh2]

    print(
        "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(
            Vh[hp.STATE].dim(), Vh[hp.PARAMETER].dim(), Vh[hp.ADJOINT].dim()
        )
    )

    extract_bottom = ExtractBottomData(mesh, Vh[hp.PARAMETER])

    #
    #  Set up the forward problem
    #
    dl.parameters["form_compiler"]["quadrature_degree"] = 3

    bc = dl.DirichletBC(Vh[hp.STATE], dl.Constant(0.0), side)

    #  Bottom and side boundary markers
    boundary_markers = dl.MeshFunction("size_t", mesh, mesh.geometry().dim() - 1)
    boundary_markers.set_all(0)

    bottom.mark(boundary_markers, 1)
    side.mark(boundary_markers, 2)
    ds = dl.Measure("ds", domain=mesh, subdomain_data=boundary_markers)

    order_ppoisson = 3.0
    functional = NonlinearPPossionForm(order_ppoisson, None, ds(1))
    pde = EnergyFunctionalPDEVariationalProblem(Vh, functional, bc, bc)

    pde.solver = dl.PETScKrylovSolver("cg", "icc")  # amg_method())
    pde.solver_fwd_inc = dl.PETScKrylovSolver("cg", "icc")  # amg_method())
    pde.solver_adj_inc = dl.PETScKrylovSolver("cg", "icc")  # amg_method())
    pde.fwd_solver.solver = dl.PETScKrylovSolver("cg", "icc")  # amg_method())

    pde.fwd_solver.parameters["gdu_tolerance"] = 1e-16
    pde.fwd_solver.parameters["LS"]["max_backtracking_iter"] = 20

    #  pde.solver.parameters["relative_tolerance"] = 1e-15
    #  pde.solver.parameters["absolute_tolerance"] = 1e-20
    #  pde.solver_fwd_inc.parameters = pde.solver.parameters
    #  pde.solver_adj_inc.parameters = pde.solver.parameters
    #  pde.fwd_solver.solver.parameters = pde.solver.parameters

    #
    # Set up the prior
    #
    gamma = 1.0
    delta = 1.0

    prior = hp.BiLaplacianPrior(Vh[hp.PARAMETER], gamma, delta, robin_bc=True)
    print(
        "Prior regularization: (delta_x - gamma*Laplacian)^order: "
        "delta={0}, gamma={1}, order={2}".format(delta, gamma, 2)
    )

    #
    #  Set up the misfit functional and generate synthetic observations
    #
    ntargets = 300
    rel_noise = 0.005

    print("Number of observation points: {0}".format(ntargets))

    if inargs["have_data"]:
        targets, data = data_file("r")
        misfit = hp.PointwiseStateObservation(Vh[hp.STATE], targets)
        misfit.d.set_local(data)

        MAX = misfit.d.norm("linf")
        noise_std_dev = rel_noise * MAX
        misfit.noise_variance = noise_std_dev * noise_std_dev
    else:
        eps = 0.05
        dummy1 = np.random.uniform(
            Length * (0.0 + eps), Length * (1.0 - eps), [ntargets, 1]
        )
        dummy2 = np.random.uniform(
            Width * (0.0 + eps), Width * (1.0 - eps), [ntargets, 1]
        )
        dummy3 = np.full((ntargets, 1), Height)
        targets = np.concatenate([dummy1, dummy2, dummy3], axis=1)

        misfit = hp.PointwiseStateObservation(Vh[hp.STATE], targets)

        mtrue = true_model(prior)

        # Export true parameter to mtrue.xdmf file
        with dl.XDMFFile(mesh.mpi_comm(), "mtrue.xdmf") as fid:
            export2XDMF(mtrue, Vh[hp.PARAMETER], fid)

        utrue = pde.generate_state()
        x = [utrue, mtrue, None]
        pde.solveFwd(x[hp.STATE], x)
        misfit.B.mult(x[hp.STATE], misfit.d)
        MAX = misfit.d.norm("linf")
        noise_std_dev = rel_noise * MAX
        misfit.noise_variance = noise_std_dev * noise_std_dev

        hp.parRandom.normal_perturb(noise_std_dev, misfit.d)

        data_file("w", target=targets, data=misfit.d.get_local())

        #  Export true state solution to uture.xdmf file
        with dl.XDMFFile(mesh.mpi_comm(), "utrue.xdmf") as fid:
            export2XDMF(utrue, Vh[hp.STATE], fid)

    model = hp.Model(pde, prior, misfit)

    #
    #  Compute the MAP point
    #
    print(sep, "Compute the MAP point", sep)

    m = prior.mean.copy()
    solver = hp.ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-8
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["max_iter"] = 25
    solver.parameters["GN_iter"] = 5
    solver.parameters["globalization"] = "LS"

    x = solver.solve([None, m, None])

    if solver.converged:
        print("\nConverged in ", solver.it, " iterations.")
    else:
        print("\nNot Converged")

    print("Termination reason:  ", solver.termination_reasons[solver.reason])
    print("Final gradient norm: ", solver.final_grad_norm)
    print("Final cost:          ", solver.final_cost)

    # Export the MAP solution to map.xdmf file
    with dl.XDMFFile(mesh.mpi_comm(), "map.xdmf") as fid:
        export2XDMF(x[hp.PARAMETER], Vh[hp.PARAMETER], fid)

    # Export the state solution corresponding to the MAP point to
    # mat_state.xdmf file
    with dl.XDMFFile(mesh.mpi_comm(), "map_state.xdmf") as fid:
        export2XDMF(x[hp.STATE], Vh[hp.STATE], fid)

    #
    #  Compute the low-rank based Laplace approximation of the posterior
    #
    print(
        sep,
        "Compute the low-rank based Laplace approximation of the posterior",
        sep,
    )
    pde.nit = 0
    model.setPointForHessianEvaluations(x, gauss_newton_approx=False)
    Hmisfit = hp.ReducedHessian(model, misfit_only=True)
    k = 100
    p = 20

    print(
        "Single/Double Pass Algorithm. Requested eigenvectors: "
        "{0}; Oversampling {1}.".format(k, p)
    )

    Omega = hp.MultiVector(x[hp.PARAMETER], k + p)
    hp.parRandom.normal(1.0, Omega)
    lmbda, V = hp.doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k)

    #plt.plot(range(0, k), lmbda, "b*", range(0, k + 1), np.ones(k + 1), "-r")
    #plt.yscale("log")
    #plt.xlabel("number")
    #plt.ylabel("eigenvalue")
    #plt.show()

    nu = hp.GaussianLRPosterior(prior, lmbda, V)
    nu.mean = x[hp.PARAMETER]

    #
    #  Set up ModPieces for implementing MCMC methods
    #
    print(sep, "Set up ModPieces for implementing MCMC methods", sep)

    # a place holder ModPiece for the parameters
    idparam = mm.IdentityOperator(Vh[hp.PARAMETER].dim())

    # log Gaussian Prior ModPiece
    gaussprior = hm.BiLaplaceGaussian(prior)
    log_gaussprior = gaussprior.AsDensity()

    # parameter to log likelihood Modpiece
    param2likelihood = hm.Param2LogLikelihood(model)

    # log target ModPiece
    log_target = mm.DensityProduct(2)

    workgraph = mm.WorkGraph()

    # Identity operator for the parameters
    workgraph.AddNode(idparam, "Identity")

    # Prior model
    workgraph.AddNode(log_gaussprior, "Prior")

    # Likelihood model
    workgraph.AddNode(param2likelihood, "Likelihood")

    # Posterior
    workgraph.AddNode(log_target, "Log_target")

    workgraph.AddEdge("Identity", 0, "Prior", 0)
    workgraph.AddEdge("Prior", 0, "Log_target", 0)

    workgraph.AddEdge("Identity", 0, "Likelihood", 0)
    workgraph.AddEdge("Likelihood", 0, "Log_target", 1)

    # Enable caching
    if inargs["MCMC"]["name"] not in ("hpcn", ""):
        log_gaussprior.EnableCache()
        param2likelihood.EnableCache()

    # Construct the problem
    postDens = workgraph.CreateModPiece("Log_target")

    mm.serveModPiece(postDens, "0.0.0.0", 4243)
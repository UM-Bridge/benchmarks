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

Input values of this script should be defined in "poisson.yaml"
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(""))))

import math
import yaml
import h5py
import numpy as np
import matplotlib.pyplot as plt


import dolfin as dl
import hippylib as hp
import muq.Modeling as mm
import muq.SamplingAlgorithms as ms
import hippylib2muq as hm


def u_boundary(x, on_boundary):
    """
    Define the boundaries that Dirichlet boundary condition is imposed on; in
    this example, they are the top and bottom boundaries.
    """
    return on_boundary and (dl.near(x[1], 0.0) or dl.near(x[1], 1.0))


class bottom_boundary(dl.SubDomain):
    """
    Define the bottom boundary.
    """

    def inside(self, x, on_boundary):
        return on_boundary and dl.near(x[1], 0.0)


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


class FluxQOI(object):
    """
    Define and evaluate the quantity of interest.
    """

    def __init__(self, Vh, dsGamma):
        self.Vh = Vh
        self.dsGamma = dsGamma
        self.n = dl.Constant((0.0, 1.0))

        self.u = None
        self.m = None
        self.L = {}

    def form(self, x):
        return (
            dl.exp(x[hp.PARAMETER])
            * dl.dot(dl.grad(x[hp.STATE]), self.n)
            * self.dsGamma
        )

    def eval(self, x):
        u = hp.vector2Function(x[hp.STATE], self.Vh[hp.STATE])
        m = hp.vector2Function(x[hp.PARAMETER], self.Vh[hp.PARAMETER])
        return np.log(dl.assemble(self.form([u, m])))


def cal_tracer(muq_samps):
    """
    Track the values of the quantity of interest.
    """
    samps_mat = muq_samps.AsMatrix()
    nums = samps_mat.shape[1]
    tracer = hp.QoiTracer(nums)

    ct = 0
    u = pde.generate_state()
    m = pde.generate_parameter()
    while ct < nums:
        m.set_local(samps_mat[:, ct])
        x = [u, m, None]
        pde.solveFwd(u, x)

        q = qoi.eval([u, m])
        tracer.append(ct, q)
        ct += 1

    return tracer


def paramcoord2eigencoord(V, B, x):
    """
    Project a parameter vector to the space spanned by eigenvectors.

    y = V^T * B * x

    :param V: eigenvectors
    :param B: the right-hand side matrix in the generalized eigenvalue problem
    :param x: parameter vector
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

if __name__ == "__main__":
    with open("poisson.yaml") as fid:
        inargs = yaml.full_load(fid)

    sep = "\n" + "#" * 80 + "\n"

    #
    # Set up the mesh and finite element function spaces
    #
    mesh = dl.UnitSquareMesh(inargs["nelement"], inargs["nelement"])

    Vh2 = dl.FunctionSpace(mesh, "Lagrange", 2)
    Vh1 = dl.FunctionSpace(mesh, "Lagrange", 1)
    Vh = [Vh2, Vh1, Vh2]

    print(
        "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(
            Vh[hp.STATE].dim(), Vh[hp.PARAMETER].dim(), Vh[hp.ADJOINT].dim()
        )
    )

    #
    # Set up the forward problem
    #
    u_bdr = dl.Expression("x[1]", degree=1)
    u_bdr0 = dl.Constant(0.0)
    bc = dl.DirichletBC(Vh[hp.STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[hp.STATE], u_bdr0, u_boundary)

    f = dl.Constant(0.0)

    def pde_varf(u, m, p):
        return (
            dl.exp(m) * dl.inner(dl.nabla_grad(u), dl.nabla_grad(p)) * dl.dx
            - f * p * dl.dx
        )

    pde = hp.PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)

    #
    # Set up the prior
    #
    gamma = 0.1
    delta = 0.5
    theta0 = 2.0
    theta1 = 0.5
    alpha = math.pi / 4
    anis_diff = dl.CompiledExpression(hp.ExpressionModule.AnisTensor2D(), degree=1)
    anis_diff.set(theta0, theta1, alpha)

    prior = hp.BiLaplacianPrior(
        Vh[hp.PARAMETER], gamma, delta, anis_diff, robin_bc=True
    )
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
        targets = np.random.uniform(0.05, 0.95, [ntargets, 2])

        misfit = hp.PointwiseStateObservation(Vh[hp.STATE], targets)

        mtrue = true_model(prior)
        if inargs["plot"]:
            objs = [
                dl.Function(Vh[hp.PARAMETER], mtrue),
                dl.Function(Vh[hp.PARAMETER], prior.mean),
            ]

            mytitles = ["True Parameter", "Prior mean"]
            hp.nb.multi1_plot(objs, mytitles)
            plt.show()

        utrue = pde.generate_state()
        x = [utrue, mtrue, None]
        pde.solveFwd(x[hp.STATE], x)
        misfit.B.mult(x[hp.STATE], misfit.d)
        MAX = misfit.d.norm("linf")
        noise_std_dev = rel_noise * MAX
        misfit.noise_variance = noise_std_dev * noise_std_dev

        hp.parRandom.normal_perturb(noise_std_dev, misfit.d)

        data_file("w", target=targets, data=misfit.d.get_local())

        if inargs["plot"]:
            vmax = max(utrue.max(), misfit.d.max())
            vmin = min(utrue.min(), misfit.d.min())
            plt.figure(figsize=(15, 5))
            hp.nb.plot(
                dl.Function(Vh[hp.STATE], utrue),
                mytitle="True State",
                subplot_loc=121,
                vmin=vmin,
                vmax=vmax,
                cmap="jet",
            )
            hp.nb.plot_pts(
                targets,
                misfit.d,
                mytitle="Observations",
                subplot_loc=122,
                vmin=vmin,
                vmax=vmax,
                cmap="jet",
            )
            plt.show()

    model = hp.Model(pde, prior, misfit)

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
    workgraph.AddNode(log_target, "Target")

    workgraph.AddEdge("Identity", 0, "Prior", 0)
    workgraph.AddEdge("Prior", 0, "Target", 0)

    workgraph.AddEdge("Identity", 0, "Likelihood", 0)
    workgraph.AddEdge("Likelihood", 0, "Target", 1)

    # Enable caching
    if inargs["MCMC"] not in ("pcn", "hpcn", ""):
        log_gaussprior.EnableCache()
        param2likelihood.EnableCache()

    # Construct the problem
    post_dens = workgraph.CreateModPiece("Target")

    mm.serveModPiece(post_dens, "posterior", "0.0.0.0", 4243)

import umbridge
import numpy as np
import cuqi


class Deconvolution1D(umbridge.Model):
    """Base benchmark for all 1D deconvolution problems"""

    dim = 128  # Dimension of the 1D signal

    def __init__(self, name):
        """Initialize the model

        parameters
        ----------
        name : str
            Name of the model

        """
        # Load data from file
        data = np.load("data_pc.npz")

        # Set up test problem
        self.TP = cuqi.testproblem.Deconvolution1D(
            dim=self.dim, phantom="pc", data=data["data"]
        )

        super().__init__(name)

    def get_input_sizes(self, config):
        return [self.dim]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        output = self.TP.posterior.logpdf(np.asarray(parameters[0]))
        return [[output]]

    def gradient(self, out_wrt, in_wrt, parameters, sens, config):
        output = self.TP.posterior.gradient(np.asarray(parameters[0]))
        return [output.tolist()]

    def supports_evaluate(self):
        return True

    def supports_gradient(self):
        return True


class Deconvolution1D_Gaussian(Deconvolution1D):
    """Deconvolution 1D with Gaussian prior"""

    def __init__(self):
        super().__init__(self.__class__.__name__)
        self.TP.prior = cuqi.distribution.GaussianCov(np.zeros(self.dim), 0.1)


class Deconvolution1D_GMRF(Deconvolution1D):
    """Deconvolution 1D with Gaussian Markov Random Field (GMRF) prior"""

    def __init__(self):
        super().__init__(self.__class__.__name__)
        self.TP.prior = cuqi.distribution.GMRF(np.zeros(self.dim), 10)


class Deconvolution1D_LMRF(Deconvolution1D):
    """Deconvolution 1D with Laplace Markov Random Field (LMRF) prior"""

    def __init__(self):
        super().__init__(self.__class__.__name__)
        self.TP.prior = cuqi.distribution.Laplace_diff(np.zeros(self.dim), 0.01)

    def supports_gradient(self):
        return False


class Deconvolution1D_CMRF(Deconvolution1D):
    """Deconvolution 1D with Cauchy Markov Random Field (CMRF) prior"""

    def __init__(self):
        super().__init__(self.__class__.__name__)
        self.TP.prior = cuqi.distribution.Cauchy_diff(np.zeros(self.dim), 0.01)


class ExactSolution(umbridge.Model):
    """Exact solution for the 1D deconvolution problem"""

    def __init__(self):
        """Initialize the model. Load the exact solution from file."""
        data = np.load("data_pc.npz")
        self.exactSolution = data["exact"]
        super().__init__(self.__class__.__name__)

    def get_input_sizes(self, config):
        return [0]

    def get_output_sizes(self, config):
        return [len(self.exactSolution)]

    def __call__(self, parameters, config):
        return [self.exactSolution.tolist()]

    def supports_evaluate(self):
        return True


model_Gaussian = Deconvolution1D_Gaussian()
model_GMRF = Deconvolution1D_GMRF()
model_LMRF = Deconvolution1D_LMRF()
model_CMRF = Deconvolution1D_CMRF()
model_exactSolution = ExactSolution()

umbridge.serve_models(
    [model_Gaussian, model_GMRF, model_LMRF, model_CMRF, model_exactSolution], 4243
)

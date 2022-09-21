import umbridge
import numpy as np
import cuqi

class Deconvolution1D(umbridge.Model):
    dim = 128

    def __init__(self):
        # Load test problem (posterior + data)
        np.random.seed(0) # Fix seed for now. Can maybe store data instead and provide it to test problem
        self.TP = cuqi.testproblem.Deconvolution1D(dim=self.dim)
        super().__init__(self.__class__.__name__)

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

model = Deconvolution1D()

umbridge.serve_models([model], 4243)

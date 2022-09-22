import umbridge
import numpy as np

# Inspired by https://github.com/chi-feng/mcmc-demo

class Donut(umbridge.Model):
    radius = 2.6
    sigma2 = 0.033

    def __init__(self):
        super().__init__("posterior")


    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        r = np.linalg.norm(parameters[0])
        return [[ - (r - Donut.radius)**2 / Donut.sigma2 ]]

    def supports_evaluate(self):
        return True

    def gradient(self, out_wrt, in_wrt, parameters, sens, config):
        r = np.linalg.norm(parameters[0])
        if (r == 0):
            return [0,0]
        return [sens[0] * parameters[0][0] * (Donut.radius / r - 1) * 2 / Donut.sigma2,
                sens[0] * parameters[0][1] * (Donut.radius / r - 1) * 2 / Donut.sigma2]

    def supports_gradient(self):
        return True

    def apply_jacobian(self, out_wrt, in_wrt, parameters, vec, config):
        r = np.linalg.norm(parameters[0])
        if (r == 0):
            return [0]
        return [vec[0] * parameters[0][0] * (Donut.radius / r - 1) * 2 / Donut.sigma2
              + vec[1] * parameters[0][1] * (Donut.radius / r - 1) * 2 / Donut.sigma2]

    def supports_apply_jacobian(self):
        return True

model = Donut()

umbridge.serve_models([model], 4243)
import umbridge
import numpy as np
from scipy.stats import multivariate_normal
import sys

# Inspired by https://github.com/chi-feng/mcmc-demo

class GaussianMixture(umbridge.Model):

    def __init__(self):
        super().__init__("posterior")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        dens1 = multivariate_normal.pdf(parameters[0], [-1.5, -1.5], 0.8)
        dens2 = multivariate_normal.pdf(parameters[0], [1.5, 1.5], 0.8)
        dens3 = multivariate_normal.pdf(parameters[0], [-2, 2], 0.5)

        if dens1 + dens2 + dens3 == 0: # log(0) not defined, so return minimal float value
            return [[ sys.float_info.min ]]
        return [[ float(np.log(dens1 + dens2 + dens3)) ]]

    def supports_evaluate(self):
        return True

    def gradient(self, out_wrt, in_wrt, parameters, sens, config):
        return [float(self.apply_jacobian(out_wrt, in_wrt, parameters, [sens[0], 0], config)[0]),
                float(self.apply_jacobian(out_wrt, in_wrt, parameters, [0, sens[0]], config)[0])]

    def supports_gradient(self):
        return True

    def apply_jacobian(self, out_wrt, in_wrt, parameters, vec, config):
        dens1 = multivariate_normal.pdf(parameters[0], [-1.5, -1.5], 0.8)
        dens2 = multivariate_normal.pdf(parameters[0], [1.5, 1.5], 0.8)
        dens3 = multivariate_normal.pdf(parameters[0], [-2, 2], 0.5)

        if dens1 + dens2 + dens3 == 0: # Return zero in log(0) case above
            return [0]

        return [float(- vec[0] / (dens1 + dens2 + dens3)
                         * (dens1 * (parameters[0][0] - -1.5) / 0.8
                         + dens2 * (parameters[0][0] - 1.5) / 0.8
                         + dens3 * (parameters[0][0] - -2) / 0.5)
                - vec[1] / (dens1 + dens2 + dens3)
                         * (dens1 * (parameters[0][1] - -1.5) / 0.8
                         + dens2 * (parameters[0][1] - 1.5) / 0.8
                         + dens3 * (parameters[0][1] - 2) / 0.5))
              ]

    def supports_apply_jacobian(self):
        return True

model = GaussianMixture()

umbridge.serve_models([model], 4243)

import umbridge
import numpy as np
from scipy.stats import multivariate_normal

# Inspired by https://github.com/chi-feng/mcmc-demo

class Banana(umbridge.Model):
    radius = 2.6
    sigma2 = 0.033

    def get_input_sizes(self):
        return [2]

    def get_output_sizes(self):
        return [1]

    def __call__(self, parameters, config={}):
        a = config.get('a', 2.0)
        b = config.get('b', 0.2)
        scale = config.get('scale', 1.0)

        y = [(parameters[0][0] / a) / scale,
             (parameters[0][1] * a + a * b * (parameters[0][0]**2 + a**2)) / scale]

        return [[multivariate_normal.logpdf(y, [0, 4], [[1.0, 0.5], [0.5, 1.0]])]]

    def supports_evaluate(self):
        return True

model = Banana()

umbridge.serve_model(model, 4243)
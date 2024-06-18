'''
https://github.com/InstituteforDiseaseModeling/EMOD-Generic/blob/main/model_covariance01/Assets/python/builder_config.py
'''

import os
import numpy as np
import umbridge
import json
from copy import deepcopy
import scipy.stats
from scipy.optimize import approx_fprime
from model import EMODForwardModel

class EMODBenchmarkModel(umbridge.Model):

    def __init__(self):
        super().__init__("posterior")
        self.model = EMODForwardModel()

    def get_input_sizes(self, config):
        return self.model.get_input_sizes(config)

    def get_output_sizes(self, config):
        return self.model.get_output_sizes(config)

    def __call__(self, parameters, config):
        # returns attack fraction of the model ([0,1]), -0 for out of bounds
        model_output = self.model(parameters, config)[0]

        model_output = [1 - model_output[0]['Susceptible Population']['Data'][-1] ]

        print(model_output)

        # enforce boundary
        if model_output[0] < 0:
            return [[-1e60]]

        # likelihood definition
        data = [0.4, 175]
        likelihood_std_dev = [0.05, 10]
        likelihood_cov_matrix_diag = [likelihood_std_dev[0]**2, likelihood_std_dev[1]**2]

        posterior = scipy.stats.multivariate_normal.logpdf(model_output, data, likelihood_cov_matrix_diag)

        return [[float(posterior)]]

    def supports_evaluate(self):
        return True    

    def gradient(self, out_wrt, in_wrt, parameters, sens, config):
        epsilon = config.get('epsilon', [0.001, 0.001, 0.001])
        grad = sens[0]*approx_fprime(parameters[0], lambda x: self([x], config)[0], epsilon=epsilon)
        return grad.tolist()

    def supports_gradient(self):
        return True            

umbridge.serve_models([EMODBenchmarkModel()], 4242)        
# umbridge.serve_models([EMODBenchmarkModel()], 4243)        

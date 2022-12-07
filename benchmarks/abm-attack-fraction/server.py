'''
https://github.com/InstituteforDiseaseModeling/EMOD-Generic/blob/main/model_covariance01/Assets/python/builder_config.py
'''

import os
import numpy as np
import umbridge
import json
from copy import deepcopy
import scipy.stats


class EMODForwardModel(umbridge.Model):
    '''
    EMOD forward model

    Todo:
        * check that running in parallel does not have a conflict.
    '''

    def __init__(self):
        super().__init__("forward")
        self.run_number = 0
        # np.random.default_rng(42)

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]


    def __call__(self, parameters, config):
        '''
        Args:
            parameters (list): [[R0], [R0 variance], [Acquisition Correlation]]
            config (dict): {}

        Returns:
            list: [[attack fraction]]
        '''
        self.run_number += 1
        model_id = deepcopy(self.run_number)
        # model_id = 1
        R0 = parameters[0][0]
        R0_VAR = parameters[0][1]
        Acquisition_Transmission_Correlation = parameters[0][2]

        config = json.load(open('config.json', 'r'))

        # check bounds
        if (R0_VAR <= 0) or (R0 <= 0) or  \
            (Acquisition_Transmission_Correlation > 1) or \
                (Acquisition_Transmission_Correlation < 0):
            # json.dump(config, open('config_{:0d}.json'.format(model_id), 'w'))
            # self._cleanup()
            return [[-np.inf]]

        inf_prd_mean   = 8.0
        inf_ln_mean    = R0/inf_prd_mean
        inf_ln_var     = R0_VAR/inf_prd_mean/inf_prd_mean
        inf_ln_sig     = np.sqrt(np.log(inf_ln_var/inf_ln_mean/inf_ln_mean+1.0))
        inf_ln_mu      = np.log(inf_ln_mean) - 0.5*inf_ln_sig*inf_ln_sig

        # load the template json config file (not thread safe)
        config['parameters']['Run_Number'] = int(model_id)
        config['parameters']['Base_Infectivity_Log_Normal_Mu'] = inf_ln_mu
        config['parameters']['Base_Infectivity_Log_Normal_Sigma'] = inf_ln_sig
        config['parameters']['Acquisition_Transmission_Correlation'] = Acquisition_Transmission_Correlation
        config['parameters']['Infectious_Period_Gaussian_Mean'] = inf_prd_mean
        json.dump(config, open('config_{:0d}.json'.format(model_id), 'w'))

        # run EMOD
        cmd = '/EMOD/Eradication -C config_{:0d}.json'.format(model_id)
        os.system(cmd)

        # read in result
        inset_chart = json.load(open(os.path.join('output', 'InsetChart.json'), 'r'))

        # attack fraction = 1 - susceptible population
        model_output = 1 - inset_chart['Channels']['Susceptible Population']['Data'][-1]

        self._cleanup()

        return [[model_output]]

    def _cleanup(self):
        '''
        Cleanup generated .json file
        '''
        # model_id = deepcopy(self.run_number)
        cmd = 'rm config_{:0d}.json'.format(self.run_number)
        os.system(cmd)

    def supports_evaluate(self):
        return True    


class EMODBenchmarkModel(umbridge.Model):
    def __init__(self):
        super().__init__("posterior")
        self.model = EMODForwardModel()

    def get_input_sizes(self, config):
        return self.model.get_input_sizes(config)

    def get_output_sizes(self, config):
        return self.model.get_output_sizes(config)

    def __call__(self, parameters, config):

        model_output = self.model(parameters, config)[0]
        # enforce boundary
        if not np.isfinite(model_output[0]):
            return [[-1e30]]

        # likelihood definition
        data = [0.4]
        likelihood_std_dev = [0.1]
        likelihood_cov_matrix_diag = [likelihood_std_dev[0]**2]

        posterior = scipy.stats.multivariate_normal.logpdf(model_output, data, likelihood_cov_matrix_diag)

        return [[posterior]]

    def supports_evaluate(self):
        return True    

umbridge.serve_models([EMODForwardModel(), EMODBenchmarkModel()], 4243)        
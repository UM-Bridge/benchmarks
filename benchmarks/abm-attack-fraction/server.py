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

class EMODForwardModel(umbridge.Model):
    '''
    EMOD forward model

    Todo:
        * check that running in parallel does not have a conflict.
    '''
    # next_id = 0

    def __init__(self):
        super().__init__("forward")
        self.run_number = 0
        # np.random.default_rng(42)
        # self.run_number = EMODForwardModel.next_id
        # EMODForwardModel.next_id += 1


    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]


    def __call__(self, parameters, config):
        '''
        Args:
            # parameters (list): [[R0], [R0 variance], [Acquisition Correlation]]
            parameters (list): [[R0], [Acquisition Correlation]]
            config (dict): {refresh_seed:True, daily_import_pressures:1.0}

        Returns:
            list: [[attack fraction]]            
        '''

        # extract config
        refresh_seed = config.get('refresh_seed', True)
        daily_import_pressures = config.get('daily_import_pressures', 1.0)
        log_level = config.get('log_level','ERROR')

        # unpack parameters
        self.run_number += 1
        model_id = deepcopy(self.run_number)
        R0 = parameters[0][0]
        R0_VAR = parameters[0][1]
        Acquisition_Transmission_Correlation = parameters[0][2]

        # load configuration and campaign files
        config_dict = json.load(open('config.json', 'r'))
        campaign_dict = json.load(open('campaign.json', 'r'))

        # check bounds
        if (R0_VAR <= 0) or (R0 <= 0) or  \
            (Acquisition_Transmission_Correlation > 1) or \
                (Acquisition_Transmission_Correlation < 0):
            return [[-999]]

        # update config file
        inf_prd_mean   = 8.0
        inf_ln_mean    = R0/inf_prd_mean
        inf_ln_var     = R0_VAR/inf_prd_mean/inf_prd_mean
        inf_ln_sig     = np.sqrt(np.log(inf_ln_var/inf_ln_mean/inf_ln_mean+1.0))
        inf_ln_mu      = np.log(inf_ln_mean) - 0.5*inf_ln_sig*inf_ln_sig
        if refresh_seed:
            config_dict['parameters']['Run_Number'] = int(model_id) % 65535
        config_dict['parameters']['Base_Infectivity_Log_Normal_Mu'] = inf_ln_mu
        config_dict['parameters']['Base_Infectivity_Log_Normal_Sigma'] = inf_ln_sig
        config_dict['parameters']['Acquisition_Transmission_Correlation'] = Acquisition_Transmission_Correlation
        config_dict['parameters']['Infectious_Period_Gaussian_Mean'] = inf_prd_mean
        config_dict['parameters'][ 'logLevel_default'] = log_level
        # update campaign file
        campaign_dict['Events'][0]['Event_Coordinator_Config']['Intervention_Config']['Daily_Import_Pressures'] = [daily_import_pressures]

        # set up campaign, config, and dempographics fils
        model_directory = f"input_{model_id:0d}"
        if not os.path.exists(model_directory):
            os.makedirs(model_directory)        
        json.dump(config_dict, open(os.path.join(model_directory, 'config.json'), 'w'))
        json.dump(campaign_dict, open(os.path.join(model_directory, 'campaign.json'), 'w'))
        cmd = f"cp demographics.json {model_directory}"
        os.system(cmd)

        # run EMOD
        cmd = f"/EMOD/Eradication -C {model_directory}/config.json -I {model_directory} -O output_{model_id:0d}"
        os.system(cmd)

        # read in result
        inset_chart = json.load(open(os.path.join(f"output_{model_id:0d}", 'InsetChart.json'), 'r'))

        # attack fraction = 1 - susceptible population
        model_output = 1 - inset_chart['Channels']['Susceptible Population']['Data'][-1]

        self._cleanup()

        return [[model_output]]

    def supports_evaluate(self):
        return True    

    def _cleanup(self):
        '''
        Cleanup generated .json file
        '''
        directory = f"input_{self.run_number:0d}"
        cmd = f"rm -rf {directory}"
        os.system(cmd)

        directory = f"output_{self.run_number:0d}"
        cmd = f"rm -rf {directory}"
        os.system(cmd)        



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

        # enforce boundary
        if model_output[0] < 0:
            return [[-1e60]]

        # likelihood definition
        data = [0.4]
        likelihood_std_dev = [0.025]
        likelihood_cov_matrix_diag = [likelihood_std_dev[0]**2]

        posterior = scipy.stats.multivariate_normal.logpdf(model_output, data, likelihood_cov_matrix_diag)

        return [[posterior]]

    def supports_evaluate(self):
        return True    

    def gradient(self, out_wrt, in_wrt, parameters, sens, config):
        epsilon = config.get('epsilon', [0.001, 0.001, 0.001])
        grad = sens[0]*approx_fprime(parameters[0], lambda x: self([x], config)[0], epsilon=epsilon)
        return grad.tolist()

    def supports_gradient(self):
        return True            

umbridge.serve_models([EMODForwardModel(), EMODBenchmarkModel()], 4243)        
# umbridge.serve_models([EMODBenchmarkModel()], 4243)        
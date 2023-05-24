'''
https://github.com/InstituteforDiseaseModeling/EMOD-Generic/blob/main/model_covariance01/Assets/python/builder_config.py
'''

import os
import numpy as np
import umbridge
import json
from copy import deepcopy

class EMODForwardModel(umbridge.Model):
    '''
    EMOD forward model

    Todo:
        * check that running in parallel does not have a conflict.
    '''

    def __init__(self):
        super().__init__("forward")
        self.run_number = 0

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]


    def __call__(self, parameters, config):
        '''
        Args:
            parameters (list): [[R0], [R0 variance], [Acquisition Correlation]]
            config (dict): {refresh_seed:True, daily_import_pressures:1.0, log_level: 'ERROR'}

        Returns:
            list: [[dict of 'New Infections', 'Infected', 'Infectious Population', 'Susceptible Population', 'Symptomatic Population', 'Recovered Population', and 'Exposed Population']]            
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

        # remove extra files
        self._cleanup()

        return [[{c:inset_chart['Channels'][c] for c in ['New Infections', 'Infected', 'Infectious Population', 'Susceptible Population', 'Symptomatic Population', 'Recovered Population', 'Exposed Population']}]]

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


if __name__ == "__main__":
    umbridge.serve_models([EMODForwardModel()], 4242)        

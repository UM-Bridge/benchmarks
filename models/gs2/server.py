import umbridge
import json
import os
import csv
from pyrokinetics import Pyro
from datetime import datetime
from fileinput import FileInput



class GS2Model(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [2] 

    def get_output_sizes(self, config):
        return [3] 

    def __call__(self, parameters, config):
        input_file = "fast.in" # Select input file
        os.system("mkdir -p restart") # GS2 needs this folder otherwise will fail
        with FileInput(files=input_file, inplace=True) as file:
            checkpoint = 0
            for line in file:
                if "&species_parameters_3" in line:
                    checkpoint += 1
                if "tprim" in line and checkpoint == 1:
                    line = f"    tprim = {parameters[0][0]}"
                if "vnewk" in line and checkpoint == 1:
                    line = f"    vnewk = {parameters[0][1]}"
                print(line.rstrip()) # FileInput redirects stdout to file
        
        # Run the model 
        mpirank = config.get("ranks", 4)
        os.system(f"mpirun -n {mpirank} /usr/gs2/bin/gs2 {input_file}")
        
        # Read results using Pyrokinetics package and print output
        gs2_input = input_file
        pyro = Pyro(gk_file=gs2_input, gk_code="GS2")
        pyro.load_gk_output()
        data = pyro.gk_output
        print(data)
        heat = data["heat"].sel(species="electron", field='phi').isel(time=-1)
        growth_rate = data["growth_rate"].isel(time=-1)
        mode_freq = data["mode_frequency"].isel(time=-1)
        output = [[heat.to_numpy().squeeze().tolist(), growth_rate.to_numpy().squeeze().tolist(), mode_freq.to_numpy().squeeze().tolist()]]
        return output 

    def supports_evaluate(self):
        return True

model = GS2Model()

umbridge.serve_models([model], 4242)

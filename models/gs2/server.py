import umbridge
import os
from pyrokinetics import Pyro



class GS2Model(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [2] 

    def get_output_sizes(self, config):
        return [3] 

    def __call__(self, parameters, config):
        input_file = "fast.in" # Select input file
        pyro = Pyro(gk_file=input_file, gk_code="GS2")
        os.system("mkdir -p restart") # GS2 needs this folder otherwise will fail
        if True: # Added to make it fast! Remove for production runs!
            pyro.gs2_input["knobs"]["nstep"] = 50
            pyro.gs2_input["theta_grid_parameters"]["ntheta"] = 10
            pyro.gs2_input["theta_grid_parameters"]["nperiod"] = 2
            pyro.gs2_input["le_grids_knobs"]["ngauss"] = 5
            pyro.gs2_input["le_grids_knobs"]["negrid"] = 2

        pyro.gs2_input["species_parameters_3"]["tprim"] = float(parameters[0][0])
        pyro.gs2_input["species_parameters_3"]["vnewk"] = float(parameters[0][1])
        pyro.gk_input.data = pyro.gs2_input
        pyro.gk_input.write(input_file)
        
        # Run the model 
        mpirank = config.get("ranks", 4)
        os.system(f"mpirun --allow-run-as-root -n {mpirank} /usr/gs2/bin/gs2 {input_file}") # --allow-run-as-root to suppress OMPI error, not reccomended when running outside of docker!
        
        # Read results and print output
        pyro = Pyro(gk_file=input_file, gk_code="GS2")
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

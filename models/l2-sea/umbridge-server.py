import umbridge
import os
import f90nml

class L2Sea(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [15]

    def get_output_sizes(self, config):
        return [5]

    def __call__(self, parameters, config):
        # Write first (and only) parameter vector to file
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'w') as f:
            for param in parameters[0][1:]:
                f.write(str(param) + '\n')
            f.close()

        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'r') as f:
            print(f.read())


        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.nml') as nml_file:
            nml = f90nml.read(nml_file)

        nml['MAIN_PARAMETERS']['igrid'] = config.get("fidelity", 7)
        nml['FREE_WARP']['fr'] = parameters[0][0]
        nml['FREE_WARP']['sinkoff'] = config.get("sinkoff", 'y')
        nml['FREE_WARP']['trimoff'] = config.get("trimoff", 'y')
        nml.write('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.aux')

        # System call, cd into working directory and call model binary
        print(config)
        os.system('cd /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415 && mv SBDF.aux SBDF.nml && ../../bin/L2-Sea')

        os.system('mv /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000/ /output')

        # Read second line of output file, split and return last 5 elements as output
        with open('/output/objective.out', 'r') as f:
            f.readline() # Skip first line
            line = f.readline() # Read second line
            line_split = line.split()[-5:] # Split and keep last 5 elements
            model_output = [float(i) for i in line_split] # Convert to float

            return [model_output]

    def supports_evaluate(self):
        return True

model = L2Sea()

umbridge.serve_models([model], 4242)

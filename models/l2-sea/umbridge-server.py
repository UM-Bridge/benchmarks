import umbridge
import os

class TestModel(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [14]

    def get_output_sizes(self, config):
        return [5]

    def __call__(self, parameters, config):
        # Write first (and only) parameter vector to file
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'w') as f:
            for param in parameters[0]:
                f.write(str(param) + '\n')

        # System call, cd into working directory and call model binary
        print(config)
        os.system('cd /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415 && ../../bin/L2-Sea')

        # Read second line of output file, split and return last 5 elements as output
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000/objective.out', 'r') as f:
            f.readline() # Skip first line
            line = f.readline() # Read second line
            line_split = line.split()[-5:] # Split and keep last 5 elements
            model_output = [float(i) for i in line_split] # Convert to float

            return [model_output]

    def supports_evaluate(self):
        return True

testmodel = TestModel()

umbridge.serve_models([testmodel], 4242)

import umbridge
import os

class TestModel(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [8]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
#    def __call__(self, parameters, config):
        arguments = " ".join([str(x) for x in parameters[0]]);

        num_threads = str(config.get("NumThreads"));
        basis_degree = str(config.get("BasisDegree"));
        fidelity = str(config.get("Fidelity"));

        arguments = arguments + " " + basis_degree + " " + fidelity

        # System call, cd into working directory and call model binary
        os.system('export IGATOOLS_NUM_THREADS=' + num_threads + ' && /build_igatools/tests/models/poisson_lorenzo/poisson_lorenzo.release ' + arguments)
#        os.system('export IGATOOLS_NUM_THREADS=10 && LD_LIBRARY_PATH=./ ./poisson_lorenzo.release ' + arguments)
#        os.system('cd ./test && export IGATOOLS_NUM_THREADS=10 && LD_LIBRARY_PATH=./ ./poisson_lorenzo.release ' + arguments)

        # Read second line of output file
        with open('/poisson_lorenzo_results.dat', 'r') as f:
#        with open('/build_igatools/tests/models/poisson_lorenzo/poisson_lorenzo_results.dat', 'r') as f:
#            f.readline() # Skip first line
            line = f.readline() # Read first line

            return [[float(line)]]


    def supports_evaluate(self):
        return True

testmodel = TestModel()

umbridge.serve_models([testmodel], 4242)

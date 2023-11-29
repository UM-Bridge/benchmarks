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
        arguments = " ".join([str(x) for x in parameters[0]]);

        num_threads = str(config.get("NumThreads"));
        basis_degree = str(config.get("BasisDegree"));
        fidelity = str(config.get("Fidelity"));

        arguments = arguments + " " + basis_degree + " " + fidelity

        # System call, cd into working directory and call model binary
        os.system('. /opt/intel/oneapi/setvars.sh && export IGATOOLS_NUM_THREADS=' + num_threads + ' && /build_igatools/tests/models/poisson_lorenzo/poisson_lorenzo.release ' + arguments)

        # Read second line of output file
        with open('/poisson_lorenzo_results.dat', 'r') as f:
            line = f.readline() # Read first line

        return [[float(line)]]


    def supports_evaluate(self):
        return True




class TestBenchmark(umbridge.Model):

    def __init__(self):
        super().__init__("benchmark")

    def get_input_sizes(self, config):
        return [8]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        arguments = " ".join([str(x) for x in parameters[0]]);

        num_threads = str(config.get("NumThreads"));
        basis_degree = str(4);
        fidelity = str(2);

        arguments = arguments + " " + basis_degree + " " + fidelity

        # System call, cd into working directory and call model binary
        os.system('. /opt/intel/oneapi/setvars.sh && export IGATOOLS_NUM_THREADS=' + num_threads + ' && /build_igatools/tests/models/poisson_lorenzo/poisson_lorenzo.release ' + arguments)

        # Read second line of output file
        with open('/poisson_lorenzo_results.dat', 'r') as f:
            line = f.readline() # Read first line

        return [[float(line)]]


    def supports_evaluate(self):
        return True





testmodel = TestModel()
testbenchmark = TestBenchmark()


umbridge.serve_models([testmodel,testbenchmark], 4242)

import umbridge
import os
import f90nml

class L2Sea_model(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [16]

    def get_output_sizes(self, config):
        return [5]

    def __call__(self, parameters, config):
        # Write first (and only) parameter vector to file
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'w') as f:
            for param in parameters[0][2:]:
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
        nml['PANCA_PARAMETERS']['ztrasla'] = parameters[0][1]

        nml.write('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.aux')

        # System call, cd into working directory and call model binary
        print(config)
        return_value = os.system('cd /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415; mv SBDF.aux SBDF.nml; ../../bin/L2-Sea')
        if return_value != 0:
            return [[0,0,0,0,0]]


        # Read second line of output file, split and return last 5 elements as output
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000/objective.out', 'r') as f:
            f.readline() # Skip first line
            line = f.readline() # Read second line
            line_split = line.split()[-5:] # Split and keep last 5 elements
            model_output = [float(i) for i in line_split] # Convert to float

            return [model_output]

    def supports_evaluate(self):
        return True

class L2Sea_benchmarkUQ(umbridge.Model):

    def __init__(self):
        super().__init__("benchmark_UQ")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        # Write first (and only) parameter vector to file
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'w') as f:
            x = [0]*14
            for i in x:
                f.write(str(i) + '\n')
            f.close()

        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'r') as f:
            print(f.read())


        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.nml') as nml_file:
            nml = f90nml.read(nml_file)

        nml['MAIN_PARAMETERS']['igrid'] = config.get("fidelity", 7)
        nml['FREE_WARP']['fr'] = parameters[0][0]
        nml['FREE_WARP']['sinkoff'] = 'y' #config.get("sinkoff", 'y')
        nml['FREE_WARP']['trimoff'] = 'y' #config.get("trimoff", 'y')
        nml['PANCA_PARAMETERS']['ztrasla'] = parameters[0][1]

        nml.write('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.aux')

        # System call, cd into working directory and call model binary
        print(config)
        return_value = os.system('cd /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415; mv SBDF.aux SBDF.nml; ../../bin/L2-Sea')
        if return_value != 0:
            return [[0]]


        # Read second line of output file, split and return last 5 elements as output
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000/objective.out', 'r') as f:
            f.readline() # Skip first line
            line = f.readline() # Read second line
            line_split = line.split()[-5:] # Split and keep last 5 elements
            model_output = [float(i) for i in line_split] # Convert to float

            return [[model_output[0]]]

    def supports_evaluate(self):
        return True

class L2Sea_benchmarkOPT(umbridge.Model):

    def __init__(self):
        super().__init__("benchmark_OPT")

    def get_input_sizes(self, config):
        return [14]

    def get_output_sizes(self, config):
        return [5]

    def __call__(self, parameters, config):
        # Write first (and only) parameter vector to file
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'w') as f:
            for param in parameters[0][0:]:
                f.write(str(param) + '\n')
            f.close()

        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/variables.inp', 'r') as f:
            print(f.read())


        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.nml') as nml_file:
            nml = f90nml.read(nml_file)

        nml['MAIN_PARAMETERS']['igrid'] = config.get("fidelity", 7)
        nml['FREE_WARP']['fr'] = 0.28 #parameters[0][0]
        nml['FREE_WARP']['sinkoff'] = 'y' #config.get("sinkoff", 'y')
        nml['FREE_WARP']['trimoff'] = 'y' #config.get("trimoff", 'y')
        nml['PANCA_PARAMETERS']['ztrasla'] = -6.16 #parameters[0][1]

        nml.write('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/SBDF.aux')

        # System call, cd into working directory and call model binary
        print(config)
        return_value = os.system('cd /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415; mv SBDF.aux SBDF.nml; ../../bin/L2-Sea')
        if return_value != 0:
            return [[0,0,0,0,0]]


        # Read second line of output file, split and return last 5 elements as output
        with open('/NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000/objective.out', 'r') as f:
            f.readline() # Skip first line
            line = f.readline() # Read second line
            line_split = line.split()[-5:] # Split and keep last 5 elements
            model_output = [float(i) for i in line_split] # Convert to float

            return [model_output]

    def supports_evaluate(self):
        return True


model = L2Sea_model()
benchmark_UQ = L2Sea_benchmarkUQ()
benchmark_OPT = L2Sea_benchmarkOPT()

umbridge.serve_models([model,benchmark_UQ,benchmark_OPT], int(os.environ.get("PORT", 4242)))

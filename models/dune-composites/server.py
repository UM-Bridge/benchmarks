import umbridge
import json
import os
import csv

class DuneCompModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [1]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):

        mpiranks = config.get("ranks",4)
        stackSeq = config.get("stack","example2.csv")

        #os.system(f"mpirun --allow-run-as-root -np {mpiranks} ./ExampleScaling -stackingSequence {stackSeq}", stdout=PIPE, stderr=PIPE)
        import subprocess
        output = subprocess.check_output(f"mpirun --allow-run-as-root --oversubscribe -np {mpiranks} ./ExampleScaling -stackingSequence {stackSeq}", shell=True)

        print(output.decode("utf-8"))

        return [[output.decode("utf-8")]]

    def supports_evaluate(self):
        return True

model = DuneCompModel()

umbridge.serve_models([model], 4242)

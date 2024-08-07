import umbridge
import json
import os
import csv
import re

class DuneCompModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [346]

    def get_output_sizes(self, config):
        return [1]

    # Input consists of KL expansion of a wrinkle
    def write_coeffs_to_file(self,filename, data):
        with open(filename, 'w') as file:
            row = data[0]
            row_str = ', '.join(f'{num:.6e}' for num in row)
            file.write(row_str + ';')

    def __call__(self, parameters, config):

        mpiranks = config.get("ranks",4)
        stackSeq = config.get("stack","example2.csv")

        self.write_coeffs_to_file("RandomFieldCoefficient/coeffs.txt", parameters)

        os.system(f"mpirun --allow-run-as-root --oversubscribe -np {mpiranks} ./ExampleScaling -stackingSequence {stackSeq} >> output.txt")
        
        with open("output.txt", 'r') as file:
            output_str = file.read()
        
        max_deflection = float(re.search(r'Maximum deflection\s*=\s*([0-9.e+-]+)', output_str).group(1))
        return [[max_deflection]]

    def supports_evaluate(self):
        return True


if __name__ == "__main__":
    umbridge.serve_models([DuneCompModel()], 4242)  


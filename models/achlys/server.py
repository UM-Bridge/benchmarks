import umbridge
import json
import os
import csv

class AchlysModel(umbridge.Model):
    def __init(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [5]

    def get_output_sizes(self, config):
        return [120]

    def __call__(self, parameters, config):

        input = {"Materials.implant.E1": parameters[0][0],
                 "Materials.implant.E2": parameters[0][1],
                 "Materials.implant.E3": parameters[0][2],
                 "Materials.implant.n1": parameters[0][3],
                 "Materials.implant.n2": parameters[0][4]}

        # Write input to JSON file
        with open("/achlys-uq/achlys-uq/input.json", "w") as f:
            f.write(json.dumps(input))

        os.system("source achlys-uq/scripts/bashrc && /achlys-uq/achlys-uq/run_desorp")

        # Read results CSV file, write rows to output
        output = []
        with open("/achlys-uq/achlys-uq/output.csv", "r") as f:
            reader = csv.reader(f)
            for row in reader:
                output.append(row)
                print(f"Read CSV Row: {row}")

        print(f"output: {output}")

        return output

    def supports_evaluate(self):
        return True

model = AchlysModel()

umbridge.serve_models([model], 4242)
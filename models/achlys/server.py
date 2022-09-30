import umbridge
import json
import os
import csv

class AchlysModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [5]

    def get_output_sizes(self, config):
        return [120]

    def __call__(self, parameters, config):

        input = f"""
            {{
                "Achlys": {{
                    "Materials": {{
                        "implant": {{
                            "E1": {parameters[0][0]},
                            "E2": {parameters[0][1]},
                            "E3": {parameters[0][2]},
                            "n1": {parameters[0][3]},
                            "n2": {parameters[0][4]}
                        }}
                    }}
                }},
                "Options": {{
                    "n_cores": 8
                }}
            }}
        """

        # Write input to JSON file
        with open("/opt/achlys-uq/achlys-uq/input.json", "w") as f:
            f.write(input)

        os.system(". /opt/achlys-uq/scripts/bashrc && cd /opt/achlys-uq/achlys-uq && /opt/achlys-uq/achlys-uq/run_desorp_umbridge")

        # Read results CSV file, write rows to output
        output = []
        with open("/opt/achlys-uq/achlys-uq/output.csv", "r") as f:
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

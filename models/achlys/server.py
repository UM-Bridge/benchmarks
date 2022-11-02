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
        return [500]

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
        with open("/opt/input.json", "w") as f:
            f.write(input)

        # add MOOSE Python libraries to PYTHONPATH, modify input files and run achlys
        os.system("cd /opt && export PYTHONPATH=$PYTHONPATH:/opt/moose/python && \
                ACHLYS_PATH=/opt/achlys/problems/thermal_desorption/ogorodnikova/tds_multiapp && \
                ./modify_input_file input.json $ACHLYS_PATH/implant_sub.i && \
                ./modify_input_file input.json $ACHLYS_PATH/resting_multi.i && \
                ./modify_input_file input.json $ACHLYS_PATH/desorp_multi.i && \
                achlys/achlys-opt --n-threads=8 -i $ACHLYS_PATH/desorp_multi.i >> /dev/null")

        # Read results CSV file, write rows to output
        output = []
        with open("/opt/achlys/problems/thermal_desorption/ogorodnikova/tds_multiapp/desorp_multi_out.csv", "r") as f:
            reader = csv.reader(f)
            for row in reader:
                output.append(row)
                print(f"Read CSV Row: {row}")

        # Interpolate data to a grid of size 500 before returning
        from scipy.interpolate import interp1d
        import pandas as pd
        from numpy import linspace, ndarray

        result = pd.read_csv('/opt/achlys/problems/thermal_desorption/ogorodnikova/tds_multiapp/desorp_multi_out.csv',usecols=['time','pfc_flux'])
        f = interp1d(result['time'],result['pfc_flux'])
        output = [ndarray.tolist(f(linspace(0,62.5,500)))]

        print(f"output: {output}")

        return output

    def supports_evaluate(self):
        return True

model = AchlysModel()

umbridge.serve_models([model], 4242)

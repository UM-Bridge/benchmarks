import umbridge
import time
import json
import os

class SolarModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [9]

    def get_output_sizes(self, config):
        return [6]

    def __call__(self, parameters, config):
        problem_id = config.get("id",1)

        with open("/solar/x.txt", "w",encoding="utf-8") as f:
            f.write('\n'.join([' '.join(str(i)) for i in parameters]))

        import subprocess
        output = subprocess.run(["/solar/bin/solar", f"""{problem_id}""", "/solar/x.txt"], capture_output=True)

        if output.returncode:
            print("Solar exited with error")
            return [[0, 0, 0, 0, 0, 0]]

        output = [float(i) for i in output.stdout.split(b" ")]
        return [output]

    def supports_evaluate(self):
        return True

model = SolarModel()

umbridge.serve_models([model], 4242)

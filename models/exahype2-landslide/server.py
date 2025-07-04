import umbridge
import time
import os


class TafjordLandslide(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [0]

    def get_output_sizes(self, config):
        return [1]
    
    def __call__(self, parameters, config):
        start = time.time()
        os.system("./TafjordLandslide.Release")
        end = time.time() - start
        return [[end]]

    def support_evaluate(self):
        return True


model = TafjordLandslide()
umbridge.serve_models([models], 4242)

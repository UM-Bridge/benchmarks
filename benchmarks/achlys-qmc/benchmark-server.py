import umbridge
import scipy.stats

class Benchmark(umbridge.Model):
    def __init__(self, model_url):
        super().__init__("posterior")
        self.model = umbridge.HTTPModel(model_url, "forward")

    def get_input_sizes(self, config):
        return self.model.get_input_sizes()

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        #TODO
        return [[posterior]]

    def supports_evaluate(self):
        return True

benchmark = Benchmark("http://localhost:4242")
forward = umbridge.HTTPModel("http://localhost:4242", "forward")

umbridge.serve_models([benchmark, forward], 4243)

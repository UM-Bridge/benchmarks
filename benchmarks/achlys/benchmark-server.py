import umbridge

import numpy
import scipy.stats as stats

from scipy.interpolate import interp1d

class Benchmark(umbridge.Model):
    def __init__(self, model_url):

        super().__init__("posterior")

        self.model = umbridge.HTTPModel(model_url, "forward")

        # set the prior distribution of parameters.
        prior_bounds = np.array([[0.7, 1.0], [0.9, 1.3], [1.1, 1.75], [5e-4, 5e-3], [1e-4, 1e-3]])
        self.prior = stats.uniform(loc=prior_bounds[:,0], scale=prior_bounds[:,1]-prior_bounds[:,0])

        # set the temperature grid of the output
        self.t = np.linspace(300, 800, 500)

        # set the data from Ogorodnikova et al. (2003)
        self.t_data = np.array([300, 310, 350, 382, 405, 430, 485, 535, 560, 575, 635, 660, 710, 720])
        self.y_data = 1e18*np.array([0.2, 0.42, 2, 3.45, 4.05, 4.1, 3.6, 2.6, 1.95, 1.5, 1, 0.55, 0.1, 0])

        # set the standard deviation and covariance of the data.
        self.sigma = 1e17
        self.cov_data = self.sigma**2 * np.eye(self.y_data.shape[0])

        # set up the likelihood function.
        self.likelihood = stats.multivariate_normal(self.y_data, self.cov_data)

    def get_input_sizes(self, config):
        return self.model.get_input_sizes()

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):

        # compute the prior log-density of the parameters.
        log_prior = self.prior.logpdf(parameters).sum()

        # evaluate the model and rescale the output from volumetric to atomic fraction.
        pfc_flux = np.array(self.model(parameters)).flatten()
        desorption_rate = 6.3e28*pfc_flux

        # create an interpolator function on the temperature grid.
        f = interp1d(self.t, desorption_rate)

        # compute the log-likelihood.
        log_like = self.likelihood.logpdf(f(self.t_data))

        # compute the posterior.
        log_posterior = log_prior + log_like

        return [[log_posterior]]

    def supports_evaluate(self):
        return True

benchmark = Benchmark("http://localhost:4242")
forward = umbridge.HTTPModel("http://localhost:4242", "forward")

umbridge.serve_models([benchmark, forward], 4243)

import umbridge
from cookiepde import CookiePDE
from dolfin import *


class CookieForward(umbridge.Model):
    """
    Class representing the forward cookie elliptic PDE model.
    """

    def __init__(self):
        """
        Initialize the CookieForward model.
        """
        super().__init__("forward")

    def get_input_sizes(self, config):
        """
        Get the sizes of the input parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the input parameter vector.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Get the sizes of the output parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the output parameter vector.
        """
        return [1]

    def __call__(self, parameters, config):
        """
        Evaluate the CookieForward model.

        Parameters:
        parameters (list): List of input parameters for the model.
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the quantity of interest.
        """
        # Fill missing entries in config with default values
        config = verifyConfig(config)
        # Initialize PDE model
        model = CookiePDE(config['N'], config['BasisDegree'])
        # Set up cookie problem
        model.setupProblem('cookie', parameters[0], config['quad_degree'], varcoeffs=config['diffzero'])
        # Solve linear system with preconditioning pc and solver tolerance tol
        model.solve(config['directsolver'], config['pc'], config['tol'])
        # Compute quantity of interest (QoI) on solution
        integral = model.computebenchmarkqoi()
        # Return QoI
        return [[integral]]

    def supports_evaluate(self):
        """
        Check if the model supports evaluation.

        Returns:
        bool: True, indicating that the model supports evaluation.
        """
        return True

    def supports_gradient(self):
        """
        Check if the model supports gradient computation.

        Returns:
        bool: False, indicating that the model does not support gradient computation.
        """
        return False


class CookieBenchmark(umbridge.Model):
    """
    Class representing the Cookie Benchmark elliptic PDE model.
    """

    def __init__(self):
        """
        Initialize the CookieBenchmark model.
        """
        super().__init__("benchmark")

    def get_input_sizes(self, config):
        """
        Get the sizes of the input parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the input parameter vector.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Get the sizes of the output parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the output parameter vector.
        """
        return [1]

    def __call__(self, parameters, config):
        """
        Evaluate the Cookie Benchmark model.

        Parameters:
        parameters (list): List of input parameters for the model.
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the quantity of interest.
        """
        # Use default values by ensuring dictionary is empty
        config={}
        config = verifyConfig(config)
        # Initialize PDE model
        model = CookiePDE(config['N'], config['BasisDegree'])
        # Set up cookie problem
        model.setupProblem(
            'cookie', parameters[0], config['quad_degree'], varcoeffs=config['diffzero'])
        # Solve linear system, possibly with preconditioning pc and solver tolerance tol
        model.solve(config['directsolver'],config['pc'], config['tol'])
        # Compute quantity of interest (QoI) on solution
        integral = model.computebenchmarkqoi()
        # Return QoI
        return [[integral]]

    def supports_evaluate(self):
        """
        Check if the model supports evaluation.

        Returns:
        bool: True, indicating that the model supports evaluation.
        """
        return True

    def supports_gradient(self):
        """
        Check if the model supports gradient computation.

        Returns:
        bool: False, indicating that the model does not support gradient computation.
        """
        return False

class CookieTime(umbridge.Model):
    """
    A model for parabolic formulation of the cookie model.

    Inherits from umbridge.Model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initializes the CookieTime object.

        Args:
            None

        Returns:
            None
        """
        super().__init__("forwardparabolic")

    def get_input_sizes(self, config):
        """
        Returns the input sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the input sizes.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Returns the output sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the output sizes.
        """
        return [1]
    
    def __call__(self, parameters, config):
        """
        Evaluates the parabolic cookie model.

        Args:
            parameters (list): A list of parameters.
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the computed integral.
        """
        # Verfiy config and fills in empty keys.
        config = verifyConfig(config)

        # Set up discrete formulation.
        model = CookiePDE(config['N'], config['BasisDegree'])
        # Note that we have an additional advection term
        model.setupProblem('cookie', parameters[0], varcoeffs=config['diffzero'], advection=config['advection'])
        # Use the custom TR-AB2 solver. Optional solveTime function uses built in PETSC TS solver.
        u = model.solveTimeSimple(config['letol'],config['T'])
        # Compute QoI at finalTime T
        integral = model.computebenchmarkqoi()
        # model.writeSln("outputFinal")
        return [[integral]]

    def supports_evaluate(self):
        """
        Checks if the model supports evaluation.

        Returns:
            bool: True if the model supports evaluation, False otherwise.
        """
        return True

    def supports_gradient(self):
        """
        Checks if the model supports gradient computation.

        Returns:
            bool: True if the model supports gradient computation, False otherwise.
        """
        return False

class CookieTimeBenchmark(umbridge.Model):
    """
    A model for parabolic formulation of the cookie model.

    Inherits from umbridge.Model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initializes the CookieTime object.

        Args:
            None

        Returns:
            None
        """
        super().__init__("benchmarkparabolic")

    def get_input_sizes(self, config):
        """
        Returns the input sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the input sizes.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Returns the output sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the output sizes.
        """
        return [1]
    
    def __call__(self, parameters, config):
        """
        Evaluates the parabolic cookie model.

        Args:
            parameters (list): A list of parameters.
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the computed integral.
        """
        # Use default config values
        config = {}
        config = verifyConfig(config)

        # Set up discrete formulation.
        model = CookiePDE(config['N'], config['BasisDegree'])
        # Note that we have an additional advection term
        model.setupProblem('cookie', parameters[0], varcoeffs=config['diffzero'], advection=config['advection'])
        # Use the custom TR-AB2 solver. Optional solveTime function uses built in PETSC TS solver.
        u = model.solveTimeSimple(config['letol'],config['T'])
        # Compute QoI at finalTime T
        integral = model.computebenchmarkqoi()
        # model.writeSln("outputFinal")
        return [[integral]]

    def supports_evaluate(self):
        """
        Checks if the model supports evaluation.

        Returns:
            bool: True if the model supports evaluation, False otherwise.
        """
        return True

    def supports_gradient(self):
        """
        Checks if the model supports gradient computation.

        Returns:
            bool: True if the model supports gradient computation, False otherwise.
        """
        return False


def verifyConfig(config):
    if config is None:
        config = {}

    # Use direct solver by default
    if 'directsolver' not in config:
        config['directsolver'] = 1

    # Use 400x400 mesh by default
    if 'Fidelity' not in config:
        config['N'] = 400
    else:
        config['N'] = int(100 * config['Fidelity'])

    # Use Q1 approximation by default
    if 'BasisDegree' not in config:
        config['BasisDegree'] = 1
    
    # Use degree 8 quadrature by default
    if 'quad_degree' not in config:
        config['quad_degree'] = 8
    
    # Use background diffusion 1.0 by default
    if 'diffzero' not in config:
        config['diffzero'] = [1.0]
    
    # Use no preconditioning by default
    if 'pc' not in config:
        config['pc'] = "none"

    # Use GM-RES tol 1e-4 by default
    if 'tol' not in config:
        config['tol'] = 1e-4
    
    # Use local timestepping error tolerance 1e-4 by default
    if 'letol' not in config:
        config['letol'] = 1e-4
    
    # Use a finalTime T = 10.0 by default
    if 'T' not in config:
        config['T'] = 10.0

    if 'advection' not in config:
        config['advection'] = 0

    print(config)
    return config

# Initialise UM-BRIDGE models
cookieforward = CookieForward()
cookiebenchmark= CookieBenchmark()
cookietimebenchmark = CookieTimeBenchmark()
cookietime = CookieTime()

# Start UM-BRIDGE server
umbridge.serve_models([cookieforward,cookiebenchmark,cookietime,cookietimebenchmark], 4242)

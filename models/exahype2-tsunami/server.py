import umbridge




class TsunamiModel(umbridge.Model):
  def __init__(self, ranks=1):
    super().__init__("forward")
    self.ranks = ranks

  def get_input_sizes(self, config):
    return [2]

  def get_output_sizes(self, config):
    return [4]

  def call(self, parameters, config):


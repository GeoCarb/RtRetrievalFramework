------------------------------------------------------------
--- Uncertainty Quantification Testing Configuration
------------------------------------------------------------

require "geocarb_base_config"

config = GeocarbBaseConfig:new()

require "uncertainty_quantification"
init_uq(config)

config:do_config()

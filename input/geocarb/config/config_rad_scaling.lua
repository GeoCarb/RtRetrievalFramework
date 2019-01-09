------------------------------------------------------------
-- Adds necessary changes for using radiance scaling
-- instead of coxmunk + lamberitan
------------------------------------------------------------

require "geocarb_base_config"

config = GeocarbBaseConfig:new()

able.insert(config.fm.instrument.instrument_correction.ic, "radiance_scaling")

config:do_config()

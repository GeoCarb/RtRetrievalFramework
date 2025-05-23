------------------------------------------------------------
--- OCO-2 is almost the same as OCO-2. Right now, the
--- differences are a different EOF file and a different spectral
--- window

require "oco_base_config"
oco3_base_config_dir = ConfigCommon.local_dir()

Oco3BaseConfig = OcoBaseConfig:new {
   static_file = oco3_base_config_dir .. "/../input/l2_oco3_static_input.h5",

   static_eof_file = oco3_base_config_dir .. "/../input/l2_oco3_eof.h5",
}

-- Use 4 EOFs for OCO3
Oco3BaseConfig.fm.instrument.instrument_correction.ic_land = { "eof_land_1", "eof_land_2","eof_land_3", "eof_land_4" }
Oco3BaseConfig.fm.instrument.instrument_correction.ic_water = { "eof_water_1", "eof_water_2","eof_water_3","eof_water_4" }

--- But only 3 EOFs for band 1 and 2
Oco3BaseConfig.fm.instrument.instrument_correction.eof_land_4.retrieve_bands = { false, false, true }
Oco3BaseConfig.fm.instrument.instrument_correction.eof_land_4.eof_used = { false, false, true }
Oco3BaseConfig.fm.instrument.instrument_correction.eof_water_4.retrieve_bands = { false, false, true }
Oco3BaseConfig.fm.instrument.instrument_correction.eof_water_4.eof_used = { false, false, true }

--- Use different absco scaling than OCO-2
Oco3BaseConfig.fm.atmosphere.absorber.CO2.table_scale = { 1.0, 0.9994, 0.9875 }
Oco3BaseConfig.fm.atmosphere.absorber.O2.table_scale = 1.001

 
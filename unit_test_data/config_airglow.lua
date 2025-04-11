------------------------------------------------------------
--- Adds the use of airglowEffect to the retrieval
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

function airglow_ap()
   local f_coeffs = Blitz_double_array_1d(2)
   f_coeffs:set(0, -1.35039e-09)
   f_coeffs:set(1, 0.0016)
   return f_coeffs
end

function airglow_cov()
   local f_cov = Blitz_double_array_2d(2,2)
   f_cov:set(Range.all(), Range.all(), 0)
   f_cov:set(0, 0, 0.02 * 0.02)
   f_cov:set(1, 1, 7e-4 * 7e-4)
   return f_cov
end

function reference_point_func()
   return DoubleWithUnit(0.755, "micron")
end

config.fm.spectrum_effect.airglow = {
   apriori = airglow_ap,
   covariance = airglow_cov,
   creator = ConfigCommon.airglow_effect,
   reference_point = reference_point_func,
   retrieved = true,
}
table.insert(config.fm.spectrum_effect.speceff, "airglow")

config:do_config()

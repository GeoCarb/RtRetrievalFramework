
--- This sets up the "standard" run we use in the unit tests, you can
--- then use this and override any feature you want to be done
--- differently
---
--- Pull anything from CommonConfig through the CommonConfig name and
--- anything in GeocarbConfig from that name instead of all through
--- GeocarbConfig to better differentiate where routines are being
--- loaded from.

require "geocarb_config"
geocarb_base_config_dir = ConfigCommon.local_dir()

GeocarbBaseConfig = GeocarbConfig:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   -- As a convenience, we pass in a few fields by environment
   -- variables.
   sid_string = os.getenv("sounding_id"),
   spectrum_file = os.getenv("spectrum_file"),
   met_file = os.getenv("met_file"),
   imap_file = os.getenv("imap_file"),
   --- Scene file is only used it we are trying to match a simulator
   --- run. So for a real data, this will be a empty string, and will
   --- in fact be ignored.
   scene_file = os.getenv("scene_file"),
   static_file = geocarb_base_config_dir .. "/../input/l2_geocarb_static_input.h5",
   static_eof_file = geocarb_base_config_dir .. "/../input/l2_geocarb_eof.h5",
   static_solar_file = config_common_dir .. "/../input/l2_solar_model.h5",
   static_aerosol_file = config_common_dir .. "/../input/l2_aerosol_combined.h5",
   -- Can have a different aerosol file for Merra aerosols
   -- static_merra_aerosol_file = config_common_dir .. "/../input/l2_aerosol_combined_RH.h5",

   merra_dir = "",

------------------------------------------------------------
-- Set this true to get diagnostic messages to help debug
-- problems with Lua
------------------------------------------------------------

   diagnostic = false,

------------------------------------------------------------
-- Paths for Absco data. We first look in local path, and if
-- we can't read the file try the full path. This allows us
-- to use a local disk on the cluster.
------------------------------------------------------------

   absco_local_path = "",
   absco_path = "",

------------------------------------------------------------
-- Connor solver
------------------------------------------------------------

   solver = { threshold=2.0,
              min_iteration=3,
              max_iteration=7,
              max_divergence=2,
              max_chisq=1.4,
              gamma_initial=10.0,
              h2o_scale_index0=20,
              h2o_scale_index1=20,
              h2o_scale_cov_initial=0.001,
              ch4_scale_index0=21,
              ch4_scale_index1=21,
              ch4_scale_cov_initial=0.001,
              co_scale_index0=22,
              co_scale_index1=22,
              co_scale_cov_initial=0.0001,
              create = ConfigCommon.connor_solver },

------------------------------------------------------------
-- Iterative solver
--
-- Thresholds and limits for the retrieval process.
-- Some thresholds are for testing the cost function.  For
-- example g_tol_abs is a tolerance for checking the
-- gradient of the cost function.  On the other hand, some
-- are used for testing the minimizer (solver).  For
-- example, minimizer_size_tol is used to check the solver.
-- Therefore, we cal this group of constants the retrieval
-- thresholds and not the solver or the problem thresholds.
------------------------------------------------------------

   -- solver =
   not_used_solver = { max_cost_function_calls=20,
                       dx_tol_abs=1e-5,
                       dx_tol_rel=1e-5,
                       g_tol_abs=1e-5,
                       minimizer_size_tol=1e-5,
                       opt_problem = ConfigCommon.nlls_max_a_posteriori,
                       iter_solver = ConfigCommon.nlls_solver_gsl_lmsder,
                       create = ConfigCommon.iterative_solver},

------------------------------------------------------------
-- If true then launch solver, otherwise just do a forward
-- model calculation, with jacobians if write_jacobians is
-- true
------------------------------------------------------------

   do_retrieval = true,

------------------------------------------------------------
-- True if we want to write the jacobians out
------------------------------------------------------------

   write_jacobian = false,

------------------------------------------------------------
-- True if we want to write high resolution spectra out
------------------------------------------------------------

   write_high_res_spectra = false,

------------------------------------------------------------
-- True if we want to generate output for every iteration
------------------------------------------------------------

   iteration_output = false,

------------------------------------------------------------
--- Log level
------------------------------------------------------------

   log_level = LogImp.INFO,

-----------------------------------------------------------
--- Default creators for everything. We default to getting
--- stuff from the HDF file, but this can be changed if
--- desired.
------------------------------------------------------------

   fm = {
      creator = ConfigCommon.oco_forward_model,
      common = {
         desc_band_name = ConfigCommon.hdf_read_string_vector("Common/desc_band_name"),
         hdf_band_name = ConfigCommon.hdf_read_string_vector("Common/hdf_band_name"),
         band_reference = ConfigCommon.hdf_read_double_with_unit_1d("Common/band_reference_point"),
         creator = ConfigCommon.table_function_eval,
      },
      spec_win = {
         creator = ConfigCommon.spectral_window_hdf,
      },
      input = {
         creator = ConfigCommon.l1b_met_input,
         l1b = {
            creator = GeocarbConfig.level1b_hdf,
            noise = {
               creator = GeocarbConfig.oco_noise,
               max_ms = { 1.4e21, 4.5e20, 2.5e20, 2.0e20 },
            },
         },
         met = {
            creator = GeocarbConfig.oco_met,
         },
      },
      stokes_coefficient = {
         creator = ConfigCommon.stokes_coefficient_constant,
         value = ConfigCommon.stokes_coefficient_l1b,
         central_wl = ConfigCommon.stokes_coefficient_central_wl_l1b,
         -- Hard code value rather than reading from l1b
--       value = ConfigCommon.stokes_coefficient_value({{{1,0}, {2,0}, {3,0}, {4,0}},
--                                                      {{5,0}, {6,0}, {7,0}, {8,0}},
--                                                      {{9,0}, {10,0},{11,0},{12,0}},
--                                                      {{13,0},{14,0},{15,0},{16,0}}}),
      },
      instrument = {
         creator = ConfigCommon.ils_instrument,
         ils_half_width = { DoubleWithUnit(4.09e-04, "um"),
                            DoubleWithUnit(1.08e-03, "um"),
                            DoubleWithUnit(1.40e-03, "um"),
                            DoubleWithUnit(1.58e-03, "um")},
         dispersion = {
            creator = ConfigCommon.dispersion_polynomial,
            apriori = ConfigCommon.l1b_spectral_coefficient_i,
            covariance = GeocarbConfig.dispersion_covariance_i("Instrument/Dispersion"),
            number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
            retrieved = true,
            is_one_based = true,
            num_parameters = 2,
         },
         ils_func = {
            creator = GeocarbConfig.ils_table_l1b,
            scale_apriori = {1.0, 1.0, 1.0, 1.0},
            scale_cov = {0.001, 0.001, 0.001, 0.001},
            retrieve_bands = {false, false, false, false},
         },
         instrument_correction = {
            creator = GeocarbConfig.instrument_correction_list_acquisition_mode,
            -- Right now we are using only one set of EOFs for all 3 modes.
            -- If we end up doing this all the time in the future, we should
            -- consider just adding a new creator that doesn't pick the EOF
            -- base on mode. But for now leave this functionality in.
            ic_nadir  = { "eof_glint_1",  "eof_glint_2",  "eof_glint_3",  "eof_glint_4"},
            ic_glint  = { "eof_glint_1",  "eof_glint_2",  "eof_glint_3",  "eof_glint_4"},
            ic_target = { "eof_glint_1",  "eof_glint_2",  "eof_glint_3",  "eof_glint_4"},
--          ic_nadir  = { "eof_nadir_1",  "eof_nadir_2",  "eof_nadir_3",  "eof_nadir_4"},
--          ic_glint  = { "eof_glint_1",  "eof_glint_2",  "eof_glint_3",  "eof_glint_4"},
--          ic_target = { "eof_target_1", "eof_target_2", "eof_target_3", "eof_target_4"},

            eof_nadir_1 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Nadir",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 1,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_nadir_2 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Nadir",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               order = 2,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_nadir_3 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Nadir",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               order = 3,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_nadir_4 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Nadir",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               order = 4,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },

            eof_glint_1 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Glint",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 1,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_glint_2 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Glint",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               order = 2,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_glint_3 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Glint",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               order = 3,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_glint_4 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Glint",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               order = 4,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },

            eof_target_1 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Target",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 1,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_target_2 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Target",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
               order = 2,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_target_3 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Target",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 3 - 1),
               order = 3,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },
            eof_target_4 = {
               hdf_group = "Instrument/EmpiricalOrthogonalFunction/Target",
               apriori = ConfigCommon.hdf_eof_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               covariance = ConfigCommon.hdf_eof_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 4 - 1),
               order = 4,
               by_pixel = true,
               scale_uncertainty = true,
               scale_to_stddev = 1e19,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = {true, true, true, true},
               eof_used = {true, true, true, true},
            },

            -- Disabled by default, add "radiance_scaling" to
            -- config.fm.instrument_correction.ic to enable.
            -- Coxmunk+Lambertian will be used instead
            radiance_scaling = {
               apriori = ConfigCommon.hdf_apriori_i("Instrument/RadianceScaling/Coxmunk"),
               covariance = ConfigCommon.hdf_covariance_i("Instrument/RadianceScaling/Coxmunk"),
               creator = ConfigCommon.radiance_scaling_sv_fit_coxmunk_only,
               retrieve_bands = {true, true, true, true},
            },
         },
      },
      spectrum_effect = {
         creator = ConfigCommon.spectrum_effect_list,
         speceff = { "solar_model", "fluorescence" },
--       speceff = { "solar_model", "instrument_doppler", "fluorescence" },
         solar_model = {
            creator = ConfigCommon.solar_absorption_and_continuum,
            doppler_shift = {
               creator = ConfigCommon.solar_doppler_from_l1b,
               do_doppler_shift = true,
            },
            solar_absorption = {
               creator = ConfigCommon.solar_absorption_table,
            },
            solar_continuum = {
               creator = ConfigCommon.solar_continuum_table,
               convert_from_photon = false,
            },
         },
         instrument_doppler = {
            creator = ConfigCommon.instrument_doppler,
            retrieved = false,
         },
         fluorescence = {
            creator = GeocarbConfig.fluorescence_effect_land_only,
            apriori = ConfigCommon.fluorescence_apriori("Fluorescence"),
            sif_sigma_scale = 1.0,
            covariance = ConfigCommon.fluorescence_covariance("Fluorescence"),
            reference_point = ConfigCommon.hdf_read_double_with_unit("Fluorescence/reference_point"),
            retrieved = true,
         },
      },
      spec_samp = {
--       creator = ConfigCommon.uniform_spectrum_sampling,
         creator = GeocarbConfig.nonuniform_spectrum_sampling,
         high_resolution_spectrum_spacing = DoubleWithUnit(0.01, "cm^-1"),
         nonunif_rt_grid_files = {
            o2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_1"),
            weak_co2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_2"),
            strong_co2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_3"),
            ch4 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_4"),
         },
      },
      rt = {
--       creator = ConfigCommon.radiative_transfer_lrad,
--       nstream = 8,
         creator = ConfigCommon.radiative_transfer_lsi,
         nadir_threshold = 1e-6,
         lsi_constant = {
            dedicated_twostream = true,
            -- Note that the "1" here is just a convention to use the
            -- dedicated two stream code
            low_stream = 1,
            -- LIDORT input is in Half-Streams. Full-streams is double
            -- this (so high_stream = 8 would mean 16 full-streams)
            high_stream = 8
         },
      },
      state_vector = {
         creator = ConfigCommon.state_vector_creator,
      },
      atmosphere = {
         creator = ConfigCommon.atmosphere_oco,
         constants = {
            creator = ConfigCommon.default_constant,
         },
         pressure = {
            apriori = ConfigCommon.met_pressure,
            covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),

--          pressure_levels = GeocarbConfig.pressure_levels_from_scene,
--          creator = ConfigCommon.pressure_fixed_level,

            a = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_a"),
            b = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_b"),
            creator = ConfigCommon.pressure_sigma,
         },
         temperature = {
            apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
            covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
            creator = ConfigCommon.temperature_met,
         },
         ground = {
            -- Instrument specific solar strengths used for ground calculations
            solar_strength = {4.876e21, 2.077e21, 1.142e21, 8.408e20},

            -- Pure lambertian
            lambertian = {
--             albedo_coeffs = {{0.2, 0.0}, {0.2, 0.0}, {0.2, 0.0}, {0.2, 0.0}},
--             albedo_coeffs = {{0.25, 0.0}, {0.30, 0.0}, {0.15, 0.0}, {0.066, 0.0}},
--             apriori = GeocarbConfig.albedo_from_constants,
               apriori = GeocarbConfig.oco_albedo_from_radiance(1),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Albedo"),
               retrieve_bands = {true, true, true, true},
               creator = ConfigCommon.lambertian_retrieval,
            },

            -- Coxmunk windspeed and refractive index inputs
            coxmunk = {
               refractive_index = ConfigCommon.hdf_apriori("Ground/Refractive_Index"),
               apriori = ConfigCommon.met_windspeed,
               covariance = ConfigCommon.hdf_covariance("Ground/Windspeed"),
               creator = ConfigCommon.coxmunk_retrieval,
            },

            -- Lambertian component of coxmunk + lambertian
            coxmunk_lambertian = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Albedo"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Albedo"),
               retrieve_bands = {true, true, true, true},
               creator = ConfigCommon.lambertian_retrieval,
            },

            -- Lambertian component of coxmunk + lambertian
            coxmunk_scaled = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Scaled"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Scaled"),
               retrieve_bands = { true, true, true, true },
               scaled_brdf_name = "CoxMunk",
               creator = ConfigCommon.brdf_scale_retrieval,
            },

            -- Brdf vegetative kernel with Rahman retrieved parameters
            brdf_veg = {
               apriori = ConfigCommon.brdf_veg_apriori("Ground/Brdf"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
               retrieve_bands = {true, true, true, true},
               creator = ConfigCommon.brdf_veg_retrieval,
            },

            -- Brdf soil kernel with Rahman retrieved parameters
            brdf_soil = {
               apriori = ConfigCommon.brdf_soil_apriori("Ground/Brdf"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
               retrieve_bands = {true, true, true, true},
               creator = ConfigCommon.brdf_soil_retrieval,
            },

--          creator = ConfigCommon.ground_lambertian,
            creator = GeocarbConfig.ground_from_ground_type,
         },
         aerosol = {
--          creator = ConfigCommon.rayleigh_only,
--          creator = ConfigCommon.aerosol_creator,
            creator = ConfigCommon.merra_aerosol_creator,
            max_aod = 0.2,
            exp_aod = 0.8,
            min_types = 2,
            max_types = 2,
            linear_aod = false,
            relative_humidity_aerosol = false,
            max_residual = 0.005,
            apriori = ConfigCommon.hdf_apriori("/Aerosol/Merra/Gaussian/Log"),
            covariance = ConfigCommon.hdf_covariance("/Aerosol/Merra/Gaussian/Log"),
            -- Lua doesn't preserve order in a table, so we have a list
            -- saying what order we want the Aerosols in
            aerosols = {"Ice", "Water", "ST"},
            Water = {
               creator = ConfigCommon.aerosol_log_shape_gaussian,
               apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
               property = ConfigCommon.hdf_aerosol_property("wc_008"),
            },
            Ice = {
               creator = ConfigCommon.aerosol_log_shape_gaussian,
               apriori_initial = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
               apriori = GeocarbConfig.tropopause_height_ap,
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
               property = ConfigCommon.hdf_aerosol_property("ice_cloud_MODIS6_deltaM_1000"),
            },
            ST = {
               creator = ConfigCommon.aerosol_log_shape_gaussian,
               apriori = function(self)
                  return ConfigCommon.lua_to_blitz_double_1d({-5.11599580975408205124,0.03,0.04})
               end,
               covariance = function(self)
                  return ConfigCommon.lua_to_blitz_double_2d({{3.24,0,0},{0,1e-8,0},{0,0,1e-4}})
               end,
               property = ConfigCommon.hdf_aerosol_property("strat"),
            },
         },
         absorber = {
            creator = ConfigCommon.absorber_creator,
            use_cache = true,
            number_sub_layers = 10,
            gases = {"CO2", "H2O", "O2", "CH4", "CO"},
            CO2 = {
--             apriori = ConfigCommon.reference_co2_apriori_met_apriori,
               apriori = ConfigCommon.tccon_co2_apriori_met,
--             apriori = GeocarbConfig.co2_apriori_from_scene,
               covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
               absco = "absco_carbo/fabiano_201803/raw/co2_carbo.hdf",
               table_scale = {1.0, 1.0, 1.0, 1.0},
               creator = ConfigCommon.vmr_level,
            },
            H2O = {
               scale_apriori = 1.0,
               scale_cov = 0.25,
               absco = "absco_carbo/fabiano_201803/raw/h2o_carbo.hdf",
               creator = ConfigCommon.vmr_met,

--             scale_apriori = 1.0,
--             scale_cov = 0.25,
--             vmr_profile = ConfigCommon.reference_h2o_apriori_met_apriori,
--             vmr_profile = GeocarbConfig.h2o_apriori_from_scene,
--             absco = "absco_carbo/fabiano_201803/raw/h2o_carbo.hdf",
--             creator = ConfigCommon.vmr_level_scaled,
            },
            O2 = {
               apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction"),
               absco = "absco_carbo/fabiano_201803/raw/o2_carbo.hdf",
               table_scale = 1.0,
               creator = ConfigCommon.vmr_level_constant_well_mixed,

--             apriori = GeocarbConfig.o2_apriori_from_scene,
--             covariance = ConfigCommon.hdf_covariance("Gas/O2"),
--             absco = "absco_carbo/fabiano_201803/raw/o2_carbo.hdf",
--             table_scale = {1.0, 1.0, 1.0, 1.0},
--             creator = ConfigCommon.vmr_level,
            },
            CH4 = {
--             apriori = ConfigCommon.reference_ch4_apriori_met_apriori,
--             apriori = GeocarbConfig.ch4_apriori_from_scene,
--             absco = "absco_carbo/fabiano_201803/raw/ch4_carbo.hdf",
--             table_scale = {1.0, 1.0, 1.0, 1.0},
--             creator = ConfigCommon.vmr_level_constant,

               scale_apriori = 1.0,
               scale_cov = 0.25,
               vmr_profile = ConfigCommon.reference_ch4_apriori_met_apriori,
--             vmr_profile = GeocarbConfig.ch4_apriori_from_scene,
               absco = "absco_carbo/fabiano_201803/raw/ch4_carbo.hdf",
               creator = ConfigCommon.vmr_level_scaled,

--             apriori = ConfigCommon.reference_ch4_apriori_met_apriori,
--             apriori = GeocarbConfig.ch4_apriori_from_scene,
--             covariance = ConfigCommon.hdf_covariance("Gas/CH4"),
--             absco = "absco_carbo/fabiano_201803/raw/ch4_carbo.hdf",
--             table_scale = {1.0, 1.0, 1.0, 1.0},
--             creator = ConfigCommon.vmr_level,
            },
            CO = {
--             apriori = ConfigCommon.reference_co_apriori_met_apriori,
--             apriori = GeocarbConfig.co_apriori_from_scene,
--             absco = "absco_carbo/fabiano_201803/raw/co_carbo.hdf",
--             table_scale = {1.0, 1.0, 1.0, 1.0},
--             creator = ConfigCommon.vmr_level_constant,

               scale_apriori = 1.0,
               scale_cov = 0.25,
               vmr_profile = ConfigCommon.reference_co_apriori_met_apriori,
--             vmr_profile = GeocarbConfig.co_apriori_from_scene,
               absco = "absco_carbo/fabiano_201803/raw/co_carbo.hdf",
               creator = ConfigCommon.vmr_level_scaled,

--             apriori = ConfigCommon.reference_co_apriori_met_apriori,
--             apriori = GeocarbConfig.co_apriori_from_scene,
--             covariance = ConfigCommon.hdf_covariance("Gas/CO"),
--             absco = "absco_carbo/fabiano_201803/raw/co_carbo.hdf",
--             table_scale = {1.0, 1.0, 1.0, 1.0},
--             creator = ConfigCommon.vmr_level,
            },
         },
         altitude = {
            number_sub_layers = 1,
            creator = ConfigCommon.hydrostatic_altitude,
         },
         relative_humidity = {
            creator = ConfigCommon.calc_relative_humidity,
         },
      },
   },
}

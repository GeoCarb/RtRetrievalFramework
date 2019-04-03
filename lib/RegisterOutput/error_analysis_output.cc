#include "error_analysis_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ErrorAnalysisOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<ErrorAnalysis>&, const blitz::Array<bool, 1>&, bool, bool, bool, bool, bool>())
REGISTER_LUA_END()
#endif

// See base class for description

void ErrorAnalysisOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  boost::function<double ()> f;
  for(int i = 0; i < err->number_spectrometer(); ++i) {
    // Skip spectrometers with missing output
    if (not spec_flag(i)) continue;

    f = boost::bind(&ErrorAnalysis::signal, err, i);
    out->register_data_source("/SpectralParameters/signal_" +
			      err->hdf_band_name(i) + "_fph", f);
    f = boost::bind(&ErrorAnalysis::noise, err, i);
    out->register_data_source("/SpectralParameters/noise_" +
			      err->hdf_band_name(i) + "_fph", f);
    f = boost::bind(&ErrorAnalysis::residual_mean_sq, err, i);
    out->register_data_source("/SpectralParameters/residual_mean_square_" +
			      err->hdf_band_name(i), f);
    f = boost::bind(&ErrorAnalysis::relative_residual_mean_sq, err, i);
    out->register_data_source("/SpectralParameters/relative_residual_mean_square_" +
			      err->hdf_band_name(i), f);
    f = boost::bind(&ErrorAnalysis::reduced_chisq, err, i);
    out->register_data_source("/SpectralParameters/reduced_chi_squared_" +
			      err->hdf_band_name(i) + "_fph", f);
  }

  out->register_data_source("/RetrievalResults/last_step_levenberg_marquardt_parameter",
			    &ErrorAnalysis::gamma_last_step, err);

  out->register_data_source("/RetrievalResults/dof_full_vector",
                            &ErrorAnalysis::degrees_of_freedom_full_vector, err);

  if (have_co2) {
      out->register_data_source("/RetrievalResults/interference_smoothing_uncert",
                                &ErrorAnalysis::interference_smoothing_uncertainty, err);
      out->register_data_source("/RetrievalResults/co2_profile_averaging_kernel_matrix",
                                &ErrorAnalysis::co2_averaging_kernel, err);
      out->register_data_source("/RetrievalResults/xco2_gain_vector",
                                &ErrorAnalysis::xco2_gain_vector, err);
      out->register_data_source("/RetrievalResults/xco2_uncert",
                                &ErrorAnalysis::xco2_uncertainty, err);
      out->register_data_source("/RetrievalResults/xco2_uncert_noise",
                                &ErrorAnalysis::xco2_uncert_noise, err);
      out->register_data_source("/RetrievalResults/xco2_uncert_smooth",
                                &ErrorAnalysis::xco2_uncert_smooth, err);
      out->register_data_source("/RetrievalResults/xco2_uncert_interf",
		                &ErrorAnalysis::xco2_uncert_interf, err);
      out->register_data_source("/RetrievalResults/dof_co2_profile",
                                &ErrorAnalysis::xco2_degrees_of_freedom, err);
      out->register_data_source("/RetrievalResults/xco2_avg_kernel",
                                &ErrorAnalysis::xco2_avg_kernel, err);
      out->register_data_source("/RetrievalResults/xco2_avg_kernel_norm",
                                &ErrorAnalysis::xco2_avg_kernel_norm, err);
      out->register_data_source("/RetrievalResults/xco2_correlation_interf",
                                &ErrorAnalysis::xco2_correlation_interf, err);
  }

  if (have_ch4) {
      out->register_data_source("/RetrievalResults/xch4_gain_vector",
                                &ErrorAnalysis::xch4_gain_vector, err);
      out->register_data_source("/RetrievalResults/xch4_uncert",
                                &ErrorAnalysis::xch4_uncertainty, err);
      out->register_data_source("/RetrievalResults/xch4_uncert_noise",
                                &ErrorAnalysis::xch4_uncert_noise, err);
      out->register_data_source("/RetrievalResults/xch4_uncert_smooth",
                                &ErrorAnalysis::xch4_uncert_smooth, err);
      out->register_data_source("/RetrievalResults/xch4_uncert_interf",
                                &ErrorAnalysis::xch4_uncert_interf, err);
      out->register_data_source("/RetrievalResults/dof_ch4_profile",
                                &ErrorAnalysis::xch4_degrees_of_freedom, err);
      out->register_data_source("/RetrievalResults/xch4_avg_kernel",
                                &ErrorAnalysis::xch4_avg_kernel, err);
      out->register_data_source("/RetrievalResults/xch4_avg_kernel_norm",
                                &ErrorAnalysis::xch4_avg_kernel_norm, err);
      if (! ch4_profile)
          ;
      else
          out->register_data_source("/RetrievalResults/ch4_profile_averaging_kernel_matrix",
                                    &ErrorAnalysis::ch4_averaging_kernel, err);
  }

  if (have_co) {
      out->register_data_source("/RetrievalResults/xco_gain_vector",
                                &ErrorAnalysis::xco_gain_vector, err);
      out->register_data_source("/RetrievalResults/xco_uncert",
                                &ErrorAnalysis::xco_uncertainty, err);
      out->register_data_source("/RetrievalResults/xco_uncert_noise",
                                &ErrorAnalysis::xco_uncert_noise, err);
      out->register_data_source("/RetrievalResults/xco_uncert_smooth",
                                &ErrorAnalysis::xco_uncert_smooth, err);
      out->register_data_source("/RetrievalResults/xco_uncert_interf",
                                &ErrorAnalysis::xco_uncert_interf, err);
      out->register_data_source("/RetrievalResults/dof_co_profile",
                                &ErrorAnalysis::xco_degrees_of_freedom, err);
      out->register_data_source("/RetrievalResults/xco_avg_kernel",
                                &ErrorAnalysis::xco_avg_kernel, err);
      out->register_data_source("/RetrievalResults/xco_avg_kernel_norm",
                                &ErrorAnalysis::xco_avg_kernel_norm, err);
      if (! co_profile)
          ;
      else
          out->register_data_source("/RetrievalResults/co_profile_averaging_kernel_matrix",
                                    &ErrorAnalysis::co_averaging_kernel, err);
  }

  out->register_data_source("/SpectralParameters/modeled_radiance",
	                    &ErrorAnalysis::modeled_radiance, err);
}

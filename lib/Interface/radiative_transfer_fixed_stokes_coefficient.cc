#include "radiative_transfer_fixed_stokes_coefficient.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

static const char* stoke_names[] = {"I", "Q", "U", "V"};
static const int max_num_stokes = 4;

// See base class for description.
Spectrum RadiativeTransferFixedStokesCoefficient::reflectance
(const SpectralDomain& Spec_domain, int Spec_index, bool Skip_jacobian) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  range_check(Spec_index, 0, number_spectrometer());
  ArrayAd<double, 2> stk;
  if(Skip_jacobian)
    stk.reference(ArrayAd<double, 2>(stokes(Spec_domain, Spec_index)));
  else
    stk.reference(stokes_and_jacobian(Spec_domain, Spec_index));

  // Send out update to Observers with spectrum index and stokes
  // values in case these need to be recorded, wraps them
  // as Spectrum for convenience and to package up any info needed
  // Only does this if any observers are registered
  if (olist.size() > 0) {
    std::vector<boost::shared_ptr<NamedSpectrum> > stoke_spectrums;
    for(int stokes_idx = 0; stokes_idx < min(stk.value().cols(), max_num_stokes); stokes_idx++) {
      Array<double, 1> stoke_val(stk.value()(Range::all(), stokes_idx));
      Array<double, 2> stoke_jac(stk.jacobian()(Range::all(), stokes_idx, Range::all()));
      SpectralRange stoke_range(ArrayAd<double, 1>(stoke_val, stoke_jac), units::inv_sr);
      boost::shared_ptr<NamedSpectrum> stoke_ptr(new NamedSpectrum(Spec_domain, stoke_range, "stokes_" + std::string(stoke_names[stokes_idx]), Spec_index));
      stoke_spectrums.push_back(stoke_ptr);
    }
    const_cast<RadiativeTransferFixedStokesCoefficient *>(this)->notify_update_do(stoke_spectrums);
  }

  ArrayAd<double, 1> res(stk.rows(), stk.number_variable());

  Array<double, 2> stokes_coef_sub(stokes_coef->stokes_coefficient().value()(Spec_index, Range(0, number_stokes() - 1), Range::all()));
/*
  res.value() = sum(stk.value()(i1, i2) * stokes_coef_sub(i2), i2);
  if(!res.is_constant()) {
    if(!stokes_coef->stokes_coefficient().is_constant()) {
      Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()(Spec_index, Range(0, number_stokes() - 1), Range::all()));
      res.jacobian() =
	sum(stk.jacobian()(i1, i3, i2) * stokes_coef_sub(i3), i3) +
	sum(stk.value()(i1, i3) * stokes_coef_jac_sub(i3, i2), i3);
    } else
      res.jacobian() =
	sum(stk.jacobian()(i1, i3, i2) * stokes_coef_sub(i3), i3);
  } else if(!stokes_coef->stokes_coefficient().is_constant()) {
      Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()(Spec_index, Range(0, number_stokes() - 1), Range::all()));
      res.jacobian() =
	sum(stk.value()(i1, i3) * stokes_coef_jac_sub(i3, i2), i3);
  }
*/
/*
  double stokes_coeff_central_wl[4] = {0., 0., 0., 0.};
*/
  double stokes_coeff_central_wl[4] = {7.650e+02, 1.606e+03, 2.065e+03, 2.323e+03};

  Array<double, 1> wl(Spec_domain.wavelength()(Range::all()));
  Array<double, 1> HV(Spec_domain.wavelength()(Range::all()));
  Array<double, 1> delta_wl(Spec_domain.wavelength()(Range::all()));
/*
  delta_wl =       wl - stokes_coeff_central_wl[Spec_index];
*/
  delta_wl = 1.0e3*wl - stokes_coeff_central_wl[Spec_index];

  ArrayAd<double, 1> myres(res.copy());

  for (int jj = 0; jj < stk.value().rows(); jj++) {
    myres.value()(jj) = stk.value()(jj,0) * (stokes_coef_sub(0,0) + stokes_coef_sub(0,1) * delta_wl(jj)) +
                        stk.value()(jj,1) * (stokes_coef_sub(1,0) + stokes_coef_sub(1,1) * delta_wl(jj)) +
                        stk.value()(jj,2) * (stokes_coef_sub(2,0) + stokes_coef_sub(2,1) * delta_wl(jj));

    myres.jacobian()(jj,Range::all()) =
                        stk.jacobian()(jj,0,Range::all()) * (stokes_coef_sub(0,0) + stokes_coef_sub(0,1) * delta_wl(jj)) +
                        stk.jacobian()(jj,1,Range::all()) * (stokes_coef_sub(1,0) + stokes_coef_sub(1,1) * delta_wl(jj)) +
                        stk.jacobian()(jj,2,Range::all()) * (stokes_coef_sub(2,0) + stokes_coef_sub(2,1) * delta_wl(jj));
  }

  return Spectrum(Spec_domain, SpectralRange(myres, units::inv_sr));
}

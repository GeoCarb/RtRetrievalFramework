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
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  range_check(Spec_index, 0, number_spectrometer());
  ArrayAd<double, 2> stk;
  if (Skip_jacobian)
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
      boost::shared_ptr<NamedSpectrum> stoke_ptr(new NamedSpectrum(Spec_domain, stoke_range,
        "stokes_" + std::string(stoke_names[stokes_idx]), Spec_index));
      stoke_spectrums.push_back(stoke_ptr);
    }
    const_cast<RadiativeTransferFixedStokesCoefficient *>(this)->notify_update_do(stoke_spectrums);
  }

  ArrayAd<double, 1> res(stk.rows(), stk.number_variable());

  Array<double, 2> stokes_coef_sub(stokes_coef->stokes_coefficient().value()
    (Spec_index, Range(0, number_stokes() - 1), Range::all()));

  Array<double, 1> delta_wl(1.0e3*Spec_domain.wavelength() -
                   stokes_coef->stokes_coefficient_central_wl().value()(Spec_index));

  Array<double, 2> factor(delta_wl.rows(), stokes_coef_sub.cols());

  int use_old = 0;
if (use_old) {
  for (int i = 0; i < factor.cols(); ++i)
    factor(Range::all(), i) = pow(delta_wl(Range::all()), (double) i);
} else {
  factor = pow(delta_wl(i1), i2);
}
if (use_old) {
  res.value() = 0.;

  for (int i = 0; i < stk.value().rows(); ++i) {
    for (int j = 0; j < stokes_coef_sub.rows(); ++j) {
      for (int k = 0; k < stokes_coef_sub.cols(); ++k)
        res.value()(i) += stk.value()(i,j) * stokes_coef_sub(j,k) * factor(i,k);
    }
  }
} else {
  res.value() = sum(sum(stk.value()(i1, i3) * stokes_coef_sub(i3, i2), i3) *
                    factor(i1, i2), i2);
}
if (use_old) {
  res.jacobian() = 0.;

  if (! res.is_constant()) {
    if (! stokes_coef->stokes_coefficient().is_constant()) {
      Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()
        (Spec_index, Range(0, number_stokes() - 1), Range::all()));

      for (int i = 0; i < stk.value().rows(); ++i) {
        for (int j = 0; j < stokes_coef_sub.rows(); ++j) {
          for (int k = 0; k < stokes_coef_sub.cols(); ++k)
            res.jacobian()(i,Range::all()) +=
              (stk.jacobian()(i,j,Range::all()) * stokes_coef_sub(j,k) +
               stk.value()(i, j) * stokes_coef_jac_sub(j,k,Range::all())) * factor(i,k);
        }
      }
    } else {
      for (int i = 0; i < stk.value().rows(); ++i) {
        for (int j = 0; j < stokes_coef_sub.rows(); ++j) {
          for (int k = 0; k < stokes_coef_sub.cols(); ++k)
            res.jacobian()(i,Range::all()) +=
              stk.jacobian()(i,j,Range::all()) * stokes_coef_sub(j,k) * factor(i,k);
        }
      }
    }
  } else if (! stokes_coef->stokes_coefficient().is_constant()) {
    Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()
      (Spec_index, Range(0, number_stokes() - 1), Range::all()));

    for (int i = 0; i < stk.value().rows(); ++i) {
      for (int j = 0; j < stokes_coef_sub.rows(); ++j) {
        for (int k = 0; k < stokes_coef_sub.cols(); ++k)
          res.jacobian()(i,Range::all()) +=
            stk.value()(i, j) * stokes_coef_jac_sub(j,k,Range::all()) * factor(i,k);
      }
    }
  }
} else {
  if (! res.is_constant()) {
    if (! stokes_coef->stokes_coefficient().is_constant()) {
      Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()
        (Spec_index, Range(0, number_stokes() - 1), Range::all()));
      res.jacobian() = sum((sum(stk.jacobian()(i1, i4, i2) * stokes_coef_sub(i4, i3), i4) +
                            sum(stk.value()(i1, i4) * stokes_coef_jac_sub(i4, i3, i2), i4)) *
                            factor(i1, i3), i3);
    } else {
      res.jacobian() = sum(sum(stk.jacobian()(i1, i4, i2) * stokes_coef_sub(i4, i3), i4) *
                           factor(i1, i3), i3);
    }
  } else if (! stokes_coef->stokes_coefficient().is_constant()) {
    Array<double, 2> stokes_coef_jac_sub(stokes_coef->stokes_coefficient().jacobian()
      (Spec_index, Range(0, number_stokes() - 1), Range::all()));
    res.jacobian() = sum(sum(stk.value()(i1, i4) * stokes_coef_jac_sub(i4, i3, i2), i4) *
                         factor(i1, i3), i3);
  }
}
  return Spectrum(Spec_domain, SpectralRange(res, units::inv_sr));
}

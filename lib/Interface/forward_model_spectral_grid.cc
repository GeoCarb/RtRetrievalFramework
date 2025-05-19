#include "forward_model_spectral_grid.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// The low resolution grid.
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::low_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectral_window->apply(inst->pixel_spectral_domain(Spec_index), Spec_index);
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// The high resolution grid, possibly nonuniform
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::high_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectrum_sampling->spectral_domain(Spec_index, low_resolution_grid(Spec_index), inst->ils_half_width(Spec_index));
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// The high resolution grid, interpolated to be uniform.
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::high_resolution_interpolated_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectrum_sampling->spectral_domain_interpolated(Spec_index, low_resolution_grid(Spec_index), inst->ils_half_width(Spec_index));
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// Pixel indexes to use for low resolution grid.
//-----------------------------------------------------------------------
  const std::vector<int> ForwardModelSpectralGrid::pixel_list(int Spec_index) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return spectral_window->grid_indexes(inst->pixel_spectral_domain(Spec_index), Spec_index);
  }

//-----------------------------------------------------------------------
/// Interpolate a spectrum to the high_resolution_interpolated_grid()
/// sampling. 
//-----------------------------------------------------------------------
Spectrum ForwardModelSpectralGrid::interpolate_spectrum
(const Spectrum& Spec_in, int Spec_index) const
{ 
  if (nus_sampling_method == 0)
    return interpolate_spectrum_linear(Spec_in, Spec_index);
  else if (nus_sampling_method == 1)
    return interpolate_spectrum_quadratic(Spec_in, Spec_index);
  else
    throw Exception("Invalid interpolation method.");
}

//-----------------------------------------------------------------------
/// Interpolate a spectrum to the high_resolution_interpolated_grid()
/// sampling. 
//-----------------------------------------------------------------------
Spectrum ForwardModelSpectralGrid::interpolate_spectrum_linear
(const Spectrum& Spec_in, int Spec_index) const
{ 
  range_check(Spec_index, 0, number_spectrometer());
  if(!spectrum_sampling->need_interpolation(Spec_index))
    return Spec_in;
  Range ra(Range::all());
  SpectralDomain hgrid_inter = high_resolution_interpolated_grid(Spec_index);
  ArrayAd<double, 1> res(hgrid_inter.data().rows(),
             Spec_in.spectral_range().data_ad().number_variable());
  Array<double, 1> spec_in_dom = Spec_in.spectral_domain().data();
  Array<double, 1> ispec_dom = 
    hgrid_inter.convert_wave(Spec_in.spectral_domain().units());
  LinearInterpolate<double, double> 
    interp(spec_in_dom.begin(),
       spec_in_dom.end(),
       Spec_in.spectral_range().data().begin(),
       LinearInterpolate<double, double>::OUT_OF_RANGE_ERROR);
  for(int i = 0; i< ispec_dom.rows(); ++i)
    res.value()(i) = interp(ispec_dom(i));
  if(res.number_variable() > 0) {
    std::vector<Array<double, 1> > yvec;
    for(int i = 0; i< spec_in_dom.rows(); ++i)
      yvec.push_back(Spec_in.spectral_range().data_ad().jacobian()(i, ra));
    LinearInterpolate<double, Array<double, 1> > 
      jacinterp(spec_in_dom.begin(), spec_in_dom.end(),
        yvec.begin(),
    LinearInterpolate<double, Array<double, 1> >::OUT_OF_RANGE_ERROR);
    for(int i = 0; i< ispec_dom.rows(); ++i)
      res.jacobian()(i,ra) = jacinterp(ispec_dom(i));
  }
  return Spectrum(hgrid_inter.clone(), 
          SpectralRange(res, Spec_in.spectral_range().units()));
}


//-----------------------------------------------------------------------
/// Quadratic interpolation.
//-----------------------------------------------------------------------
Array<double, 1> interpolate_quadratic(const Array<double, 1>& x, const Array<double, 1>& y, const Array<double, 1>& xx, bool bound)
{
    int i = 0;
    int p = 0;
    int nx  = x.rows();
    int nxx = xx.rows();

    Array<double, 1> yy(nxx);

    for ( ; i < nxx; ++i) {
        double this_xx = xx(i);

        while (x(p) < this_xx && p < nx - 1)
            p++;
        p--;
        if (p < 0) p = 0;

        if (this_xx == x(p))
            yy(i) = y(p);
        else {
            if (p < 1) {
                yy(i) = y(p) + (y(p+1) - y(p)) * (this_xx - x(p)) / (x(p+1) - x(p));
/*
                double f_1 = (this_xx - x(p+1)) / (x(p)   - x(p+1));
                double f_2 = (this_xx - x(p))   / (x(p+1) - x(p));

                yy(i) = f_1 * y(p) + f_2 * y(p+1);
*/
            }
            else {
                double x_1  = x(p-1);
                double x_2  = x(p);
                double x_3  = x(p+1);
                double xx_i = this_xx;

                double f_1 = (xx_i - x_2) * (xx_i - x_3) / ((x_1 - x_2) * (x_1 - x_3));
                double f_2 = (xx_i - x_1) * (xx_i - x_3) / ((x_2 - x_1) * (x_2 - x_3));
                double f_3 = (xx_i - x_1) * (xx_i - x_2) / ((x_3 - x_1) * (x_3 - x_2));

                yy(i) = f_1 * y(p-1) + f_2 * y(p) + f_3 * y(p+1);
            }
        }
    }

    if (bound) {
        for (int i = 0; i < nxx; ++i) {
            if (xx(i) < x(0))      yy(i) = y(0);
            if (xx(i) > x(nx - 1)) yy(i) = y(nx - 1);
        }
    }

    return yy;
}


//-----------------------------------------------------------------------
/// Interpolate a spectrum to the high_resolution_interpolated_grid()
/// sampling. 
//-----------------------------------------------------------------------
Spectrum ForwardModelSpectralGrid::interpolate_spectrum_quadratic
(const Spectrum& Spec_in, int Spec_index) const
{
  range_check(Spec_index, 0, number_spectrometer());
  if(!spectrum_sampling->need_interpolation(Spec_index))
    return Spec_in;

  Range ra(Range::all());

  SpectralDomain hgrid_inter = high_resolution_interpolated_grid(Spec_index);
  ArrayAd<double, 1> res(hgrid_inter.data().rows(),
             Spec_in.spectral_range().data_ad().number_variable());
  Array<double, 1> spec_in_dom = Spec_in.spectral_domain().data();
  Array<double, 1> ispec_dom = 
    hgrid_inter.convert_wave(Spec_in.spectral_domain().units());

  res.value()(ra) = interpolate_quadratic(spec_in_dom, Spec_in.spectral_range().data(), ispec_dom, true);

  if(res.number_variable() > 0) {
    std::vector<Array<double, 1> > yvec;
    for(int i = 0; i< spec_in_dom.rows(); ++i)
      yvec.push_back(Spec_in.spectral_range().data_ad().jacobian()(i, ra));

    for (int i = 0; i < Spec_in.spectral_range().data_ad().number_variable(); ++i)
        res.jacobian()(ra, i) = interpolate_quadratic(spec_in_dom, Spec_in.spectral_range().data_ad().jacobian()(ra, i), ispec_dom, true);
  }

  return Spectrum(hgrid_inter.clone(), 
          SpectralRange(res, Spec_in.spectral_range().units()));
}

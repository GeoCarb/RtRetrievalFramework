#ifndef OCO_SIM_APRIORI_H
#define OCO_SIM_APRIORI_H
#include "pressure.h"
#include "hdf_file.h"
#include "oco_sounding_id.h"
#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to calculate a CO2 apriori using the scene
  file from a OCO simulator run.
*******************************************************************/

class OcoSimApriori : public Printable<OcoSimApriori> {
public:
  OcoSimApriori(const std::string& Oco_sim_scene,
		const HdfSoundingId& Sid);
  virtual ~OcoSimApriori() {}

  double h2o_vmr(double P) const { return interp_h2o(P); }
  double o2_vmr(double P) const { return interp_o2(P); }
  double co2_vmr(double P) const { return interp_co2(P); }
  double ch4_vmr(double P) const { return interp_ch4(P); }
  double co_vmr(double P) const { return interp_co(P); }

  blitz::Array<double, 1> pressure_levels() const;
  blitz::Array<double, 1> h2o_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> o2_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> co2_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> ch4_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> co_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> surface_albedo(int i) const;

  void print(std::ostream& Os) const { Os << "OcoSimApriori"; }

private:
  blitz::Array<double, 1> plev;

  LinearInterpolate<double, double> interp_h2o;
  LinearInterpolate<double, double> interp_o2;
  LinearInterpolate<double, double> interp_co2;
  LinearInterpolate<double, double> interp_ch4;
  LinearInterpolate<double, double> interp_co;
};
}

#endif

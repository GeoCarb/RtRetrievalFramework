// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "oco_sim_apriori.h"
%}
%base_import(generic_object)
%import "hdf_sounding_id.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::OcoSimApriori);
namespace FullPhysics {
class OcoSimApriori : public GenericObject {
public:
  OcoSimApriori(const std::string& Oco_sim_scene,
		const HdfSoundingId& Sid);
  double co2_vmr(double P) const;
  blitz::Array<double, 1> pressure_levels() const;
  blitz::Array<double, 1> h2o_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> o2_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> co2_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> ch4_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> co_vmr_grid(const Pressure& P) const;
  blitz::Array<double, 1> surface_albedo(int i) const;
  std::string print_to_string() const;
};
}

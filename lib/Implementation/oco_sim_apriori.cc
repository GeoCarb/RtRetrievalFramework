#include "oco_sim_apriori.h"
#include "hdf_file.h"
#include "oco_sounding_id.h"
#include <boost/algorithm/string.hpp>

using namespace FullPhysics;
using namespace blitz;
using namespace boost::algorithm;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(OcoSimApriori)
.def(luabind::constructor<const std::string&,
     const HdfSoundingId&>())
.def("pressure_levels", &OcoSimApriori::pressure_levels)
.def("h2o_vmr_grid", &OcoSimApriori::h2o_vmr_grid)
.def("o2_vmr_grid", &OcoSimApriori::o2_vmr_grid)
.def("co2_vmr_grid", &OcoSimApriori::co2_vmr_grid)
.def("ch4_vmr_grid", &OcoSimApriori::ch4_vmr_grid)
.def("co_vmr_grid", &OcoSimApriori::co_vmr_grid)
.def("surface_albedo", &OcoSimApriori::surface_albedo)
REGISTER_LUA_END()
#endif

OcoSimApriori::OcoSimApriori(const std::string& Oco_sim_scene,
			     const HdfSoundingId& Sid)
{
  HdfFile h(Oco_sim_scene);
  // We take a HdfSoundingId because that works better with Lua,
  // but we require this to actually be a OcoSoundingId.
  int row = (dynamic_cast<const OcoSoundingId&>(Sid)).scene_index();
  int nlay = h.read_field<int, 1>("/Simulation/Thermodynamic/num_layers",
				  TinyVector<int, 1>(row),
				  TinyVector<int, 1>(1))(0);

  plev.resize(nlay + 1);
  plev = h.read_field<double, 2>("/Simulation/Thermodynamic/pressure_level",
			         TinyVector<int, 2>(row, 0),
			         TinyVector<int, 2>(1, nlay + 1))(0, Range::all());

  TinyVector<int, 2> sz = h.read_shape<2>("/Simulation/Gas/species_id");
  Array<std::string, 1> species_id(h.read_field<std::string, 2>
				   ("/Simulation/Gas/species_id",
				    TinyVector<int, 2>(row, 0),
				    TinyVector<int, 2>(1, sz[1]))
				   (0, Range::all()));
  int air_dry_i = -1;
  int h2o_i = -1;
  int o2_i  = -1;
  int co2_i = -1;
  int ch4_i = -1;
  int co_i = -1;
  for(int i = 0; i < species_id.rows(); ++i) {
    std::string t = species_id(i);
    trim(t);
    if(t == "AIR_dry")
      air_dry_i = i;
    else if(t == "H2O")
      h2o_i = i;
    else if(t == "O2")
      o2_i = i;
    else if(t == "CO2")
      co2_i = i;
    else if(t == "CH4")
      ch4_i = i;
    else if(t == "CO")
      co_i = i;
  }
  if(air_dry_i == -1)
    throw Exception("AIR_dry not found in file.");
  if(h2o_i == -1)
    throw Exception("H2O not found in file.");
  if(o2_i == -1)
    throw Exception("O2 not found in file.");
  if(co2_i == -1)
    throw Exception("CO2 not found in file.");
  if(ch4_i == -1)
    throw Exception("CH4 not found in file.");
  if(co_i == -1)
    throw Exception("CO not found in file.");

  TinyVector<int, 3> sz2 = h.read_shape<3>("/Simulation/Gas/species_density");
  Array<double, 2> sden = h.read_field<double, 3>
    ("/Simulation/Gas/species_density",
     TinyVector<int, 3>(row, 0, 0),
     TinyVector<int, 3>(1, sz2[1], nlay))(0, Range::all(), Range::all());

  Array<double, 1> air_dry(sden(air_dry_i, Range::all()));

  Array<double, 1> pvmr(air_dry.shape());
  for(int i = 0; i < pvmr.rows(); ++i)
    pvmr(i) = (plev(i) + plev(i + 1)) / 2;

  Array<double, 1> vmr(air_dry.shape());

  Array<double, 1> h2o(sden(h2o_i, Range::all()));
  vmr = h2o / air_dry;
  interp_h2o = LinearInterpolate<double, double>(pvmr.begin(), pvmr.end(),
                                                 vmr.begin());

  Array<double, 1> o2(sden(o2_i, Range::all()));
  vmr = o2 / air_dry;
  interp_o2 = LinearInterpolate<double, double>(pvmr.begin(), pvmr.end(),
                                                 vmr.begin());

  Array<double, 1> co2(sden(co2_i, Range::all()));
  vmr = co2 / air_dry;
  interp_co2 = LinearInterpolate<double, double>(pvmr.begin(), pvmr.end(),
                                                 vmr.begin());

  Array<double, 1> ch4(sden(ch4_i, Range::all()));
  vmr = ch4 / air_dry;
  interp_ch4 = LinearInterpolate<double, double>(pvmr.begin(), pvmr.end(),
                                                 vmr.begin());
  Array<double, 1> co(sden(co_i, Range::all()));
  vmr = co / air_dry;
  interp_co = LinearInterpolate<double, double>(pvmr.begin(), pvmr.end(),
                                                 vmr.begin());
}

//-----------------------------------------------------------------------
/// Return CO2 VMR for each pressure in the given pressure grid.
//-----------------------------------------------------------------------

blitz::Array<double, 1> OcoSimApriori::pressure_levels() const
{
  Array<double, 1> pressure_levels(plev);
  return pressure_levels;
}

blitz::Array<double, 1> OcoSimApriori::h2o_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = h2o_vmr(p(i));
  return res;
}

blitz::Array<double, 1> OcoSimApriori::o2_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = o2_vmr(p(i));
  return res;
}

blitz::Array<double, 1> OcoSimApriori::co2_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = co2_vmr(p(i));
  return res;
}

blitz::Array<double, 1> OcoSimApriori::ch4_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = ch4_vmr(p(i));
  return res;
}

blitz::Array<double, 1> OcoSimApriori::co_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = co_vmr(p(i));
  return res;
}

blitz::Array<double, 1> OcoSimApriori::surface_albedo(int i) const
{
  Array<double, 1> albedo(2);

  cerr << "ERROR: OcoSimApriori::surface_albedo() not implemented yet\n";
  exit(1);

  return albedo;
}

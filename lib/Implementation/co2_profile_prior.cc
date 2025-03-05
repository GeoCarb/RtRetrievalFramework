#include "co2_profile_prior.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(CO2ProfilePrior)
.def(luabind::constructor<const OcoMetFile&,
     const HdfFile&>())
.def(luabind::constructor<const OcoMetFile&,
     const HdfFile&, const std::string&>())
.def(luabind::constructor<const AcosMetFile&,
     const HdfFile&>())
.def(luabind::constructor<const AcosMetFile&,
     const HdfFile&, const std::string&>())
.def("apriori_vmr", &CO2ProfilePrior::apriori_vmr)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
CO2ProfilePrior::CO2ProfilePrior
(const OcoMetFile& Met_file,
 const HdfFile& Profile_file)
: model_press(Met_file.pressure_levels())
{
  boost::shared_ptr<HdfSoundingId> hsid = Met_file.sounding_id();
  std::string field = "CO2Prior/co2_prior_profile_cpr";
  TinyVector<int, 3> sz = Profile_file.read_shape<3>(field);
  Array<double, 3> traw = Profile_file.read_field<double, 3>
    (field,
     TinyVector<int, 3>(hsid->frame_number(), hsid->sounding_number(), 0),
     TinyVector<int, 3>(1,1,sz[2]));
  co2_vmr.resize(traw.extent(thirdDim));
  co2_vmr = traw(0, 0, Range::all());
}

CO2ProfilePrior::CO2ProfilePrior
(const OcoMetFile& Met_file,
 const HdfFile& Profile_file,
 const std::string& field)
: model_press(Met_file.pressure_levels())
{
  boost::shared_ptr<HdfSoundingId> hsid = Met_file.sounding_id();
  TinyVector<int, 3> sz = Profile_file.read_shape<3>(field);
  Array<double, 3> traw = Profile_file.read_field<double, 3>
    (field,
     TinyVector<int, 3>(hsid->frame_number(), hsid->sounding_number(), 0),
     TinyVector<int, 3>(1,1,sz[2]));
  co2_vmr.resize(traw.extent(thirdDim));
  co2_vmr = traw(0, 0, Range::all());
}

CO2ProfilePrior::CO2ProfilePrior
(const AcosMetFile& Met_file,
 const HdfFile& Profile_file)
: model_press(Met_file.pressure_levels())
{
  boost::shared_ptr<HdfSoundingId> hsid = Met_file.sounding_id();
  std::string field = "CO2Prior/co2_prior_profile_cpr";
  TinyVector<int, 4> sz = Profile_file.read_shape<4>(field);
  Array<double, 4> traw = Profile_file.read_field<double, 4>
    (field,
     TinyVector<int, 4>(hsid->frame_number(), hsid->sounding_number(), 0, 0),
     TinyVector<int, 4>(1,1,1,sz[3]));
  co2_vmr.resize(traw.extent(fourthDim));
  co2_vmr = traw(0, 0, 0, Range::all());
}

CO2ProfilePrior::CO2ProfilePrior
(const AcosMetFile& Met_file,
 const HdfFile& Profile_file,
 const std::string& field)
: model_press(Met_file.pressure_levels())
{
  boost::shared_ptr<HdfSoundingId> hsid = Met_file.sounding_id();
  TinyVector<int, 4> sz = Profile_file.read_shape<4>(field);
  Array<double, 4> traw = Profile_file.read_field<double, 4>
    (field,
     TinyVector<int, 4>(hsid->frame_number(), hsid->sounding_number(), 0, 0),
     TinyVector<int, 4>(1,1,1,sz[3]));
  co2_vmr.resize(traw.extent(fourthDim));
  co2_vmr = traw(0, 0, 0, Range::all());
}

//-----------------------------------------------------------------------
/// Apriori value
//-----------------------------------------------------------------------

blitz::Array<double, 1> CO2ProfilePrior::apriori_vmr
(const Pressure& pressure) const
{
  Array<double, 1>
    press_levels(pressure.pressure_grid().convert(units::Pa).value.value());
  LinearInterpolate<double, double>
    mod_vmr_interp(model_press.begin(), model_press.end(), co2_vmr.begin());

  Array<double, 1> interp_vmr(press_levels.shape());
  for(int i = 0; i < interp_vmr.rows(); i++)
    interp_vmr(i) = mod_vmr_interp(press_levels(i));
  return interp_vmr;
}


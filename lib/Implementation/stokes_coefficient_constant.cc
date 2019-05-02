#include "stokes_coefficient_constant.h"
#include "ostream_pad.h"

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(StokesCoefficientConstant, StokesCoefficient)
.def(luabind::constructor<const blitz::Array<double, 3>&>())
.def(luabind::constructor<const blitz::Array<double, 3>&,
                          const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
StokesCoefficientConstant::StokesCoefficientConstant
(const blitz::Array<double, 3>& Stokes_coeff)
{
  stokes_coeff.resize(Stokes_coeff.shape(), 0);
  stokes_coeff = Stokes_coeff;
  stokes_coeff_central_wl.resize(Stokes_coeff.shape()(0), 0);
  stokes_coeff_central_wl = 0.;
  // This is a dummy value, since the base SubStateVectorArray expects
  // at least one value. But we mark this as unused, and don't do
  // anything with these values.
  blitz::Array<double, 1> coeff(1);
  blitz::Array<bool, 1> used(1);
  coeff(0) = 0.0;
  used(0) = false;
  init(coeff, used);
}

StokesCoefficientConstant::StokesCoefficientConstant
(const blitz::Array<double, 3>& Stokes_coeff,
 const blitz::Array<double, 1>& Stokes_coeff_central_wl)
{
  stokes_coeff.resize(Stokes_coeff.shape(), 0);
  stokes_coeff = Stokes_coeff;
  stokes_coeff_central_wl.resize(Stokes_coeff.shape()(0), 0);
  stokes_coeff_central_wl = Stokes_coeff_central_wl;
  // This is a dummy value, since the base SubStateVectorArray expects
  // at least one value. But we mark this as unused, and don't do
  // anything with these values.
  blitz::Array<double, 1> coeff(1);
  blitz::Array<bool, 1> used(1);
  coeff(0) = 0.0;
  used(0) = false;
  init(coeff, used);
}

void StokesCoefficientConstant::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "StokesCoefficientConstant:\n";
  opad << stokes_coefficient().value();
  opad.strict_sync();
}

boost::shared_ptr<StokesCoefficient> StokesCoefficientConstant::clone() const
{
  return boost::shared_ptr<StokesCoefficient>
    (new StokesCoefficientConstant(stokes_coefficient().value()));
}

#ifndef ERROR_ANALYSIS_OUTPUT_H
#define ERROR_ANALYSIS_OUTPUT_H
#include "register_output_base.h"
#include "error_analysis.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ErrorAnalysis class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ErrorAnalysis class.
*******************************************************************/
class ErrorAnalysisOutput : public RegisterOutputBase {
public:
  ErrorAnalysisOutput(const boost::shared_ptr<ErrorAnalysis>& E, 
		      const blitz::Array<bool, 1>& Spec_flag,
		      bool Have_co2 = false,
		      bool ch4_profile = false,
		      bool co_profile = false)
    : err(E), spec_flag(Spec_flag), have_co2(Have_co2), ch4_profile(ch4_profile),
      co_profile(co_profile) {}

  virtual ~ErrorAnalysisOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  mutable boost::shared_ptr<ErrorAnalysis> err;
  bool have_co2;
  bool ch4_profile;
  bool co_profile;
  const blitz::Array<bool, 1> spec_flag;
};
}
#endif

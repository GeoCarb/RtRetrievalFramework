// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This file was auto-generated

%{
#include "lidort_interface_masters.h"
%}

namespace FullPhysics {



class Brdf_Lin_Sup_Masters {

public:
  Brdf_Lin_Sup_Masters();
  
  Brdf_Sup_Inputs& brdf_sup_in();
  Brdf_Linsup_Inputs& brdf_linsup_in();
  Brdf_Input_Exception_Handling& brdf_sup_inputstatus();
  Brdf_Sup_Outputs& brdf_sup_out();
  Brdf_Linsup_Outputs& brdf_linsup_out();
  Brdf_Output_Exception_Handling& brdf_sup_outputstatus();
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Brdf_Sup_Masters {

public:
  Brdf_Sup_Masters();
  
  Brdf_Sup_Inputs& brdf_sup_in();
  Brdf_Input_Exception_Handling& brdf_sup_inputstatus();
  Brdf_Sup_Outputs& brdf_sup_out();
  Brdf_Output_Exception_Handling& brdf_sup_outputstatus();
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Lidort_Inputs {

public:
  Lidort_Inputs();
  
  Lidort_Sup_Inout& lidort_sup();
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Input_Exception_Handling& lidort_inputstatus();
  
  void sup_init();
  void read_config(const std::string& filnam_in);
  void sleave_sup_init();
  void ss_sup_init();
  void brdf_sup_init();
};


class Lidort_L_Inputs {

public:
  Lidort_L_Inputs();
  
  Lidort_Linsup_Inout& lidort_linsup();
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Fixed_Lininputs& lidort_linfixin();
  Lidort_Modified_Lininputs& lidort_linmodin();
  Lidort_Input_Exception_Handling& lidort_inputstatus();
  
  void linsup_init();
  void ss_linsup_init();
  void sleave_linsup_init();
  void brdf_linsup_init();
  void read_config(const std::string& filnam_in);
};


class Lidort_Lcs_Masters {

public:
  Lidort_Lcs_Masters();
  
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Sup_Inout& lidort_sup();
  Lidort_Outputs& lidort_out();
  Lidort_Fixed_Lininputs& lidort_linfixin();
  Lidort_Modified_Lininputs& lidort_linmodin();
  Lidort_Linsup_Inout& lidort_linsup();
  Lidort_Linoutputs& lidort_linout();
  
  void run();
};


class Lidort_Lps_Masters {

public:
  Lidort_Lps_Masters();
  
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Sup_Inout& lidort_sup();
  Lidort_Outputs& lidort_out();
  Lidort_Fixed_Lininputs& lidort_linfixin();
  Lidort_Modified_Lininputs& lidort_linmodin();
  Lidort_Linsup_Inout& lidort_linsup();
  Lidort_Linoutputs& lidort_linout();
  
  void run();
};


class Lidort_Masters {

public:
  Lidort_Masters();
  
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Sup_Inout& lidort_sup();
  Lidort_Outputs& lidort_out();
  
  void run();
};


class Lidort_Brdf_Sup_Accessories {

public:
  Lidort_Brdf_Sup_Accessories(boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_in, boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_in);
  
  Brdf_Sup_Outputs& brdf_sup_out();
  Lidort_Fixed_Inputs& lidort_fixin();
  Lidort_Modified_Inputs& lidort_modin();
  Lidort_Sup_Inout& lidort_sup();
  
  void set_brdf_inputs();
};

}
# TNANP_project

To run the TNANP_project.py must simply first select which tasks one wants to do by selecting the boolean "True" at the beggining of the script:

TASK: do_non_relativistic_polytropic_EoS                      
Description: Computes pressure and mass functions of a non-relativistic neutron star described by a polytropic equation of state.
Returns: Plot of Figure 8 from "Irina Sagert et al 2006 Eur. J. Phys. 27 577"

TASK: do_all_regimes_EoS      
Description: Computes pressure and mass functions of a neutron star with an equation of state valid for all regimes.
Returns: Plot of Figures 9 and 10 from "Irina Sagert et al 2006 Eur. J. Phys. 27 577"

TASK: do_all_regimes_EoS_specific_p0 '\n'
Description: Computes pressure and mass functions of a neutron star with an equation of state valid for all regimesc with a specific dimensionless initial defined by the variable "P0". The default choice is 0.01.
Returns: Plot of Figure 5 from "R. R. Silbar, S. Reddy American Journal of Physics 72, 892 (2004)"

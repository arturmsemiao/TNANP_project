import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, LogLocator, AutoLocator
from TNANP_module import get_star_mass_radius
import argparse

#Tasks to be done: 
do_non_relativistic_polytropic_EoS = True
do_all_regimes_EoS = True
do_all_regimes_EoS_specific_p0 = True

####################################################################################
################# NON RELATIVISTIC POLYTROPIC EQUATION OF STATE ####################
####################################################################################


if do_non_relativistic_polytropic_EoS: #In this section I intend to treat a polytropic non relativistic equation of state

    central_pressures = np.logspace(31.2,33.8,50) #dyne/cm**2
    
    _,newton_masses, newton_final_radii = get_star_mass_radius(
        central_pressures, polytropic = True, relativistic = False, TOV = False
    )
    _,TOV_masses, TOV_final_radii = get_star_mass_radius(
        central_pressures, polytropic = True, relativistic = False, TOV = True
    )
    
    #The "get_star_mass_radius" returns, for a given set of central pressure, all 
    #values of the mass and pressure functions. To get only the final masses, one
    #must select the last value of the array for each central pressure  which is
    #the correct full mass. 
    
    newton_final_masses = np.array([mass[-1] for mass in newton_masses.values()]) 
    TOV_final_masses = np.array([mass[-1] for mass in TOV_masses.values()]) 

    #FINAL RADII AND MASSES VS CENTRAL PRESSURES PLOT
    fig, ax1 = plt.subplots(figsize = (10,8))
    ax1.set_title('Mass and Radius as a function of initial central pressures $p_0$ ')
    ax1.plot(central_pressures, newton_final_radii.values(), label = 'R Newton', color = 'tab:orange')
    ax1.plot(central_pressures, TOV_final_radii.values(), linestyle='dashed',  label = 'R TOV', color = 'tab:orange')
    ax1.set_ylim(15,33)
    ax1.set_ylabel('R in km')
    ax1.legend(loc = 'upper center')
    ax2 = ax1.twinx()
    ax2.plot(central_pressures, newton_final_masses, label = 'M Newton', color = 'blue')
    ax2.plot(central_pressures, TOV_final_masses, linestyle='dashed', label = 'M TOV', color = 'blue')
    ax2.set_ylim(0,0.55)
    ax2.set_ylabel(u"Mass in $M_\u2609$")
    ax2.legend(loc = 'upper right')
    ax1.yaxis.set_minor_locator(MultipleLocator(1))
    ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax2.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.set_xlim(2e31,4e33)                                                             
    ax1.set_xlabel('$p_0$ in dyne/$cm^2$')
    ax1.set_xscale('log')
    plt.savefig('NeutronStar_polytropic_EoS_radius_vs_pressure.png')
    plt.close()    

###################################################################################################
###############################  ALL REGIMES EQUATION OF STATE  ###################################
###################################################################################################

if do_all_regimes_EoS: #In this section I intend to treat an equation of state valid for all relativistic regimes

    central_pressures = np.logspace(30,41,100) # dyne/cm**2
    
    _,newton_masses, newton_final_radii = get_star_mass_radius(
        central_pressures, polytropic = False, relativistic = False, TOV = False
    )
    _,TOV_masses, TOV_final_radii = get_star_mass_radius(
        central_pressures, polytropic = False, relativistic = False, TOV = True
    )

    #The "get_star_mass_radius" returns, for a given set of central pressure, all 
    #values of the mass and pressure functions. To get only the final masses, one
    #must select the last value of the array for each central pressure  which is
    #the correct full mass. 
    
    newton_final_masses = np.array([mass[-1] for mass in newton_masses.values()]) 
    TOV_final_masses = np.array([mass[-1] for mass in TOV_masses.values()]) 

    #FINAL RADII AND MASSES VS CENTRAL PRESSURES PLOT
    fig, ax1 = plt.subplots(figsize = (10,8))
    ax1.set_title('Mass and Radius as a function of initial central pressures $p_0$ ')
    ax1.plot(central_pressures, newton_final_radii.values(), label = 'R Newton', color = 'tab:orange')
    ax1.plot(central_pressures, TOV_final_radii.values(), linestyle='dashed',  label = 'R TOV', color = 'tab:orange')
    ax1.set_ylim(0,47)
    ax1.set_ylabel('R in km')
    ax1.legend(loc = 'upper center')
    ax2 = ax1.twinx()
    ax2.plot(central_pressures, newton_final_masses, label = 'M Newton', color = 'blue')
    ax2.plot(central_pressures, TOV_final_masses, linestyle='dashed', label = 'M TOV', color = 'blue')
    ax2.set_ylim(0,0.8)
    ax2.set_ylabel(u"Mass in $M_\u2609$")
    ax2.legend(loc = 'upper right')
    ax1.yaxis.set_minor_locator(MultipleLocator(5))
    ax1.yaxis.set_major_locator(MultipleLocator(10))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax2.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.set_xlim(1e30,1e41)                                                             
    ax1.set_xlabel('$p_0$ in dyne/$cm^2$')
    ax1.set_xscale('log')
    plt.savefig('NeutronStar_mass_radius_vs_pressure.png')
    plt.close()    

    # FINAL RADII VS FINAL MASSES PARAMETRIC PLOT
    plt.figure()
    plt.plot(TOV_final_radii.values(),TOV_final_masses)
    plt.xlim(0,26)
    plt.xlabel('R in km')
    plt.ylim(0,1)
    plt.ylabel('Mass in $M_\u2609$')
    plt.savefig('mass_vs_radius.png')
    plt.close()


#############################################################################################################
############################ FULL RELATIVISTIC (FR) TOV vs NEWTON p0 = 0.01 #################################
#############################################################################################################

if do_all_regimes_EoS_specific_p0: #In this section I intend to treat the same equation of state but now with a specifica initial pressure

    KAPPA = 5.48846e36 # dyne/cm**2 or ergs/cm**3, kappa factor for the numerically fitted EoS.

    dimless_central_pressures = np.array([0.01]) #
    central_pressures = dimless_central_pressures*KAPPA # Render dimension by multiplying by kappa factor
    
    newton_dimless_pressures, newton_masses, newton_final_radii = get_star_mass_radius(
        central_pressures, polytropic = False, relativistic = False, TOV = False
    )
    TOV_dimless_pressures, TOV_masses, TOV_final_radii = get_star_mass_radius(
        central_pressures, polytropic = False, relativistic = False, TOV = True
    )

    #Now we only have one key in the dictionary that correponds to the dimless central pressure p0 = 0.01

    #Build radius arrays to plot mass and pressure vs radius to study how these behave:
    newton_radius_array = np.linspace(0, (len(newton_masses[0.01])-1)/100, len(newton_masses[0.01]))
    TOV_radius_array = np.linspace(0, (len(TOV_masses[0.01])-1)/100, len(TOV_masses[0.01]))
    
    plt.figure()
    plt.plot(newton_radius_array, newton_dimless_pressures[0.01], label = f'$p_0$ = 0.01 Newton')
    plt.plot(TOV_radius_array, TOV_dimless_pressures[0.01], label = f'$p_0$ = 0.01 TOV')
    plt.xlabel('Radius in Km')
    plt.xlim(0,15)
    plt.ylabel('Dimless Pressure')
    plt.legend()
    plt.savefig('pressure_vs_radius.png')
    plt.close()

    plt.figure()
    plt.plot(newton_radius_array, newton_masses[0.01], label = f'$p_0$ = 0.01 Newton')
    plt.plot(TOV_radius_array, TOV_masses[0.01], label = f'$p_0$ = 0.01 TOV')
    plt.ylabel('M in $M_\u2609$')
    plt.xlim(0,15)
    plt.ylim(0,1.1)
    plt.xlabel('R in km')
    plt.legend()
    plt.savefig('mass_vs_radius.png')
    plt.close()

    print('Newton radius with p0 = 0.01: ', newton_final_radii[0.01])
    print('Newton mass with p0 = 0.01: ', newton_masses[0.01][-1])
    print('TOV radius with p0 = 0.01: ', TOV_final_radii[0.01])
    print('TOV mass with p0 = 0.01: ', TOV_masses[0.01][-1])


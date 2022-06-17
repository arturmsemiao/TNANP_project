import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, LogLocator, AutoLocator
import numpy as np
import math

R_0 = 1.476 #Half of the Schwarzchild radius in km
DYNE_CM2_TO_MSUNC2_KM3 = 5.623e-40 #Convertion factor from dyne/cm**2 to Msun*c**2/km**3 units

#Definition of constants used in equations of state
GAMMA_NONREL, GAMMA_REL = 5/3, 1 
K_NONREL, K_REL = 6.432e-26, 1/3
ALPHA_NONREL, ALPHA_REL = 1, 3*R_0
A_NONRELA, A_RELA = 2.4216, 2.8863

def newton_structure_equation(radius, dimless_energy_density, dimless_mass):
    """
    Computes newton's structure function for the pressure of the star
    """
    return - R_0 * dimless_mass/(radius**2) * dimless_energy_density

def tov_structure_equation(radius, dimless_energy_density, dimless_mass, dimless_pressure, beta):
    """
    Computes TOV's structure function for the pressure of the star
    """
    first_correction = dimless_mass/radius**2 + beta*radius*dimless_pressure 
    second_correction = 1+dimless_pressure/dimless_energy_density
    third_correction = (1-2*R_0*dimless_mass/radius)
    
    return -R_0*dimless_energy_density*first_correction*second_correction/third_correction

def mass_structure_equation(radius, dimless_energy_density, beta):
    """
    Computes structure function for the mass of the star
    """
    return beta * radius**2 * dimless_energy_density

def coupled_structure_equations(radius, dimless_pressure_and_mass, beta, equation_of_state, TOV):
    """
    Returns the output of the coupled system of differential equations for the
    pressure and mass of the star at a specific radius given an equation of state
    and a boolean TOV condition.

    dimless_pressure_and_mass is an array of length 2 where the pressure is the 1st
    element and the mass is the 2nd
    """
    dimless_pressure, dimless_mass = dimless_pressure_and_mass
    
    if dimless_pressure > 0:

        dimless_energy_density = equation_of_state(dimless_pressure)
        
        if TOV: dp_dr = tov_structure_equation(radius, dimless_energy_density, dimless_mass, dimless_pressure, beta)
        else: dp_dr = newton_structure_equation(radius, dimless_energy_density, dimless_mass)
    
        dm_dr = mass_structure_equation(radius, dimless_energy_density, beta)

    else: dp_dr, dm_dr = 0,0

    return np.array([dp_dr, dm_dr])

def runge_kutta_4(h, radius, dimless_pressure, dimless_mass, beta, equation_of_state, TOV):
    """
    Returns the new values for the pressure and mass functions of the star after
    a step h by computing the factors of the Runge Kutta 4th order model evaluated
    on the coupled system of structure equations
    """
    dimless_pressure_and_mass = np.array([dimless_pressure, dimless_mass])
    
    k1 = coupled_structure_equations(radius, dimless_pressure_and_mass, beta, equation_of_state, TOV)
    k2 = coupled_structure_equations(radius + h/2, dimless_pressure_and_mass + (h/2)*k1, beta, equation_of_state, TOV)
    k3 = coupled_structure_equations(radius + h/2, dimless_pressure_and_mass + (h/2)*k2, beta, equation_of_state, TOV)
    k4 = coupled_structure_equations(radius + h, dimless_pressure_and_mass + h*k3, beta, equation_of_state, TOV)

    new_dimless_pressure, new_dimless_mass = dimless_pressure_and_mass + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    return new_dimless_pressure, new_dimless_mass

def get_star_mass_radius(central_pressures, polytropic, relativistic, TOV):
    """
    Returns: three dictionaries with the values being the pressure function, the mass function and
    the final radius of the neutron star for each central pressure which are the keys of the dictionaries 

    Arguments:

    central_pressures: Array of initial central pressures (in dyne/cm**2) for which one wants to
    compute the neutron star's mass and radius

    polytropic: boolean. If True, one uses polytropic EoS in structure equations, if False, one uses
    the numerical general fit.
    
    relativistic: boolean. If True, one uses the relativistic polytropic EoS, if False, one uses the
    non relativistic polytropic EoS.

    TOV: boolean. If True, one uses GR corrected structure equations to compute neutrons star's mass
    and radius, if False, one neglects GR corrections and uses newtonian.
    """
    #Create dictionares that will associate to each central pressure the corresponding mass, pressure and radius arrays of the star 
    central_pressure_to_dimless_pressure_values,  central_pressure_to_final_star_radii, central_pressure_to_dimless_mass_values = {}, {}, {}
    
    MAX_RK4_STEPS = 100000 
    RK4_STEP_SIZE = 0.01

    INITIAL_RADIUS = 1e-5
    INITIAL_MASS = 0

    #Choose equation of state (EoS) between the polytropic (relativistic or non relativistic) or the general fit
    if polytropic:
        
        if relativistic: polytropic_factor, gamma, alpha = K_REL, GAMMA_REL, ALPHA_REL 
        else: polytropic_factor, gamma, alpha = K_NONREL, GAMMA_NONREL, ALPHA_NONREL

        kappa = ((1/polytropic_factor) * (R_0/alpha)**(gamma) )**(1/(gamma-1)) #ergs/cm3 or dyne/cm2
        
        def polytropic_EoS(dimless_pressure):
            return (dimless_pressure/(polytropic_factor*kappa**(gamma-1)))**(1/gamma)

        equation_of_state = polytropic_EoS
            
    else:

        kappa = 5.48846e36 #ergs/cm3 or dyne/cm2
        
        def general_fit_EoS(dimless_pressure):
            return A_NONRELA*dimless_pressure**(3/5) + A_RELA*dimless_pressure
    
        equation_of_state = general_fit_EoS
        
    beta = 4*math.pi*kappa*DYNE_CM2_TO_MSUNC2_KM3 #1/km**3, beta factor that is present in structure equations
    dimless_central_pressures = central_pressures/kappa #Render central pressures dimensionless by dividing by the kappa factor
    
    for central_pressure in dimless_central_pressures: #iterate over various initial dimensionless central pressures p0

        #Set the initial conditions of our problem
        dimless_pressure, dimless_mass, star_radius = central_pressure, INITIAL_MASS, INITIAL_RADIUS
        #Create arrays with initial conditions to register all values of pressure and mass for each step of the star's radius
        dimless_pressure_values, dimless_mass_values = np.array([dimless_pressure]), np.array([dimless_mass])
        
        step = 0

        while dimless_pressure > 0.0 and step <= MAX_RK4_STEPS:            
            
            #Update pressure and mass values with Runge Kutta 4 method
            dimless_pressure, dimless_mass = runge_kutta_4(RK4_STEP_SIZE, star_radius, dimless_pressure, dimless_mass, beta, equation_of_state, TOV)
            #Update radius with step size of RK4 model
            star_radius += RK4_STEP_SIZE 
            
            #Append to respective arrays
            dimless_pressure_values = np.append(dimless_pressure_values, dimless_pressure)
            dimless_mass_values = np.append(dimless_mass_values, dimless_mass)

            step += 1

        #Cut out last value since it corresponds to an already negative pressure
        dimless_pressure_values = dimless_pressure_values[0:-1]
        dimless_mass_values = dimless_mass_values[0:-1]

        #Registar star's final radius by subtracting the final RK4 step which was added already when pressure < 0
        final_star_radius = round(star_radius - RK4_STEP_SIZE, 4)
            
        #Record pressure and mass functions and final radius of the star with given initial dimensionless central pressure
        central_pressure_to_dimless_pressure_values[central_pressure] = dimless_pressure_values
        central_pressure_to_dimless_mass_values[central_pressure] = dimless_mass_values
        central_pressure_to_final_star_radii[central_pressure] = final_star_radius
        
    return central_pressure_to_dimless_pressure_values, central_pressure_to_dimless_mass_values, central_pressure_to_final_star_radii

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import numpy as np
from PIL import Image
import numpy as np
import rebound
import streamlit.components.v1 as components
    
class orbital_evolution():
    """ Class that has 5 methods which contribute to tracking the orbital evolution of a satellite with respect to gravitational perturbations from closeby neighbors (albeit in the same constellation, now this version of ORBYTES)
    
    Attributes
    __________
    file: str
        Pandas DataFrame that holds all sets of orbital element data
    orbital elements: floats
        a, e, inc, asc, per, true_anomaly
     
    Methods
    _______
    extract_orbital_elements()
        Extract orbital elements (ex. a, e, inc, arg of perigee, etc.) for each satellite from respective pandas DataFrame
    eccentric_anomaly()
        Calculate the eccentric anomaly by iterating over combinations of the mean_anomale & eccentricity of individual satellite
    plot()
        Plot 3 graphs which help visualize the deviation in semi-major axis (a), eccentricity (e), and inclination (inc) that a satellite undergoes in the presence of gravitational perturbers from a small cohort of neighboring bodies
    simulate_system()
        Track how the orbital elements (a, e, and inc) for an indexed body in a system of len(data) particles change in magnitude of the course of N orbits. Extract a (Noutputs x 4) dimension table of data.
    orbital_evolution_sat()
        Evolve system of particles while keeping track of the evolution of a, e, and inc for one indexed particle. Plot the results to aid in our visualization of the element's fluctuations
    """
    def __init__(self):
        return

    def eccentric_anomaly(self, file, index = 0, iterations = 500):
        file = pd.read_csv(file)
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        
        if mean_anomaly >= 0.6:
            initial_guess = np.pi
        else:
            initial_guess = mean_anomaly
        ind = 0
        while ind < iterations:
            E_old = initial_guess
            E_new = E_old - ((E_old - eccentricity*np.sin(E_old)-mean_anomaly) / (1-eccentricity*np.cos(E_old)))
            
            if np.abs(E_new - E_old) <= 1e-16:
                return E_new 
            else:
                pass
    
            initial_guess = E_new
            ind+=1
            
    def extract_orbital_elements(self, file, index = 0):
        og_file = file
        file = pd.read_csv(file)
        # Calls 'plot_orbit' function using the TLE elements defined in output file
        row_orbital_elements = file.loc[file.index[index]]
        object_name = file['OBJECT_NAME'].iloc[index]
        semi_major_axis = file['SEMIMAJOR_AXIS'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        inclination = file['INCLINATION'].iloc[index]
        right_ascension = file['RA_OF_ASC_NODE'].iloc[index]
        argument_pericenter = file['ARG_OF_PERICENTER'].iloc[index] # Angle between lines of nodes (ascending node) and periapsis 
        
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        ecc_anomaly = self.eccentric_anomaly(og_file)
        true_anomaly = np.arctan((np.sqrt(1 - eccentricity**2) * np.sin(ecc_anomaly)) / (np.cos(ecc_anomaly) - eccentricity))# Solve for true anomaly using mean anomaly and eccentricity : defines the position of a body along the ellipse at a specific time ('epoch') - angle between periapsis & current location along orbit
    
        return object_name, semi_major_axis, eccentricity, inclination, right_ascension, argument_pericenter, mean_anomaly
        
    def plot(self, a,b,c,d, object_name, orbits):

        rows, cols = 3, 1
        fig, (ax1, ax2, ax3) = plt.subplots(rows, cols, sharex = True)
        ax1.plot((((a/(2*np.pi))*92.8)/60)/24, c)
        ax1.set_ylabel('Semi-major Axis', fontsize = 8)
        ax1.tick_params(labelsize = 8)
        ax1.set_xlim(0,orbits)
        
        ax2.plot((((a/(2*np.pi))*92.8)/60)/24, b)
        ax2.set_ylabel('Eccentricity', fontsize = 8)
        ax2.tick_params(labelsize = 7.5)
        ax2.set_xlim(0,orbits)
        
        ax3.plot((((a/(2*np.pi))*92.8)/60)/24, d * (360/(2*np.pi)))
        ax3.set_ylabel('Inclination', fontsize = 8)
        ax3.tick_params(labelsize = 8)
        ax3.set_xlabel('Time (Days)')
        ax3.set_xlim(0, (orbits*92.8)/(60*24))

        ax1.set_title(f'Orbital Evolution of Satellite {object_name}', fontsize = 10, pad = 20)
        plt.plot()

    def simulate_system(self, file, index, orbits, Noutputs):
        data = pd.read_csv(file)
        
        sim = rebound.Simulation()
        sim.add('Earth') # Primary particle
        for i in range(len(data)):
            object_name, semi_major_axis, eccentricity, inclination, right_ascension, argument_pericenter, mean_anomaly = self.extract_orbital_elements(file, index = i)
            mean_anomaly = np.radians(mean_anomaly)
            f = rebound.M_to_f(eccentricity, mean_anomaly)
            sim.add(a = semi_major_axis / 6790 , e = eccentricity, inc = np.radians(inclination), Omega = np.radians(right_ascension), omega = np.radians(argument_pericenter), f = f)
        particles = sim.particles
        sim.integrator = 'ias15'
        
        torb = 2 * np.pi # One period ~ 1 Orbit around Earth ~ 92.8 minutes
        times = np.linspace(0, orbits * torb, Noutputs)
        sim.dt = 1e-3
        
        a = np.zeros(Noutputs)
        e = np.zeros(Noutputs)
        inc = np.zeros(Noutputs)
        
        for i, time in enumerate(times):
            sim.integrate(time, exact_finish_time = 0)
            orbit = sim.particles[index].calculate_orbit(primary = sim.particles[0])
            a[i] = orbit.a
            e[i] = orbit.e
            inc[i] = orbit.inc
        
        orbital_prop_SS = np.vstack((times, e, a, inc))
        orbital_evolution = np.transpose(orbital_prop_SS)
        
        return orbital_evolution
    
    def orbital_evolution_sat(self, data, index, orbits, Noutputs):
        orbital_evolution = self.simulate_system(data, index, orbits, Noutputs)
        df = pd.DataFrame(orbital_evolution, columns = ['times(orbits)', 'e', 'a', 'inc'])
        object_name_final, a, e, i, asc, per, mean = self.extract_orbital_elements(data, index)
        self.plot(df['times(orbits)'], df['e'], df['a'], df['inc'], object_name =object_name_final, orbits = orbits)
        

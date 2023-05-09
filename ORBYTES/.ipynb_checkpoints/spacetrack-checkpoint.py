import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import numpy as np
from PIL import Image
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import mpld3
import streamlit.components.v1 as components
import plotly.express as px
%matplotlib qt

class Satellite_Modeler():
    
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
        
        
    def extract_orbital_elements_sat(self, file, index):
        whole_file = file
        file = pd.read_csv(file)
        # Extracts orbital elements within each two-line element set (TLE) that is uniquely defined for satellite
        row_orbital_elements = file.loc[file.index[index]]
        object_name = file['OBJECT_NAME'].iloc[index]
        semi_major_axis = file['SEMIMAJOR_AXIS'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        inclination = file['INCLINATION'].iloc[index]
        right_ascension = file['RA_OF_ASC_NODE'].iloc[index]
        argument_pericenter = file['ARG_OF_PERICENTER'].iloc[index] 
        
        # All calulated (or extracted from data) at stated epoch
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        ecc_anomaly = self.eccentric_anomaly(whole_file)
        # Solve for true anomaly using mean anomaly and eccentricity 
        # : defines the position of a body along the ellipse at a specific time ('epoch') 
        # - angle between periapsis & current location along orbit
        true_anomaly = np.arctan((np.sqrt(1 - eccentricity**2) * np.sin(ecc_anomaly)) / (np.cos(ecc_anomaly) - eccentricity))
        
        return object_name, semi_major_axis, eccentricity, inclination, right_ascension, argument_pericenter, true_anomaly
        
        
    def earth(self, earth_radius = 6378 , e_earth = 0): # in km
        # Create a respresentation using vectors of Earth's model size
        unit_vector = np.array([1,1,1])
        earth_pos_vector = []
        for cartesian_component in unit_vector:
            earth_pos_vector.append(earth_radius / cartesian_component)
        r_x = earth_pos_vector[0]
        r_y = earth_pos_vector[1]
        r_z = earth_pos_vector[2]
        
        # Define spherical properties of Earth
        phi = np.linspace(0, np.pi, 100) # in radians
        theta = np.linspace(0, 2*np.pi, 100)
        
        # Transform to spherical coordinates 
        x = r_x * np.outer(np.cos(theta), np.sin(phi))
        y = r_y * np.outer(np.sin(theta), np.sin(phi))
        z = r_z * np.outer(np.ones_like(theta), np.cos(phi))
        
        # Find extent of Earth's orbital plane as reference plane (equatorial plane)
        r_earth = (earth_radius * (1 - e_earth**2)) / (1 + e_earth * np.cos(theta))
        
        # Polar coordinates of Earth's Equatorial Plane projection
        polar_x_e = r_earth * np.cos(theta)
        polar_y_e = r_earth * np.sin(theta)
        polar_z_e = r_earth * 0
        
        position_earth_3d = np.matrix(list(zip(polar_x_e, polar_y_e, polar_z_e)))
        
        return r_earth, x, y, z, polar_x_e, polar_y_e, polar_z_e
        
    
    def orientation_sat(self, inc = 0, asc = 0, per = 0):
        # Transformation b/w equatorial coordinates and ecliptic coordinates == 
        # most solar system body's orbits & apparent positions lie in the ecliptic plane

        # Must define rotational matrices that will properly orient the ellipse from equatorial to ecliptic coordinate system
        inc = inc * ((2*np.pi) / 360)  # in radians
        Minc = np.matrix([[1, 0, 0], # rotation about the x-axis (actually, line of nodes)
                      [0, np.cos(inc), -np.sin(inc)], 
                      [0, np.sin(inc), np.cos(inc)]])
        asc = asc * ((2*np.pi) / 360) # in radians
        Masc = np.matrix([[np.cos(asc), -np.sin(asc), 0],
                        [np.sin(asc), np.cos(asc), 0],
                        [0, 0, 1]]) # rotation about the Z-axis of plane of reference 
        per = per * ((2*np.pi) / 360) # in radians
        Mper = np.matrix([[np.cos(per), -np.sin(per), 0], 
                        [np.sin(per), np.cos(per), 0],
                       [0, 0, 1]]) # Rotation about the Z-axis of plane of orbit (actually, angular momentum vector)
        
        return Minc, Masc, Mper # Set of Rotation Matrices that orient satellite's trajectory
    
    def polar_eqn_ellipse(self, a, e):
        theta = np.linspace(0, 2* np.pi, 100)
        # Kepler's First Law -> orbits are ellipses -> DEFINE polar equation of ellipse (foundation of orbit structure)
        # We are returned radius of ellipse (which is the orbit of satellite) + polar coordinates
        r = (a * (1 - e**2)) / (1 + e * np.cos(theta))
        polar_x = r * np.cos(theta)
        polar_y = r * np.sin(theta)
        polar_z = (0) * theta
        
        return r, polar_x, polar_y, polar_z
    
    def orbit_pos(self, a, e, inc, asc, per, true_anomaly, satellite_body = False):
        # Take satellite's ORBITAL ELEMENTS & reconstruct elliptical orbit in 3D
        theta = np.linspace(0, 2 * np.pi, 100)
        Minc, Masc, Mper = self.orientation_sat(inc, asc, per)
        r, polar_x, polar_y, polar_z = self.polar_eqn_ellipse(a, e)
        
        # Co-planar satellite trajectory with Earth's equator (has yet to be oriented appropriately)
        if satellite_body == False:
            pos_sat_orbit = np.matrix(list(zip(polar_x, polar_y, polar_z)))
        else:
            satellite_angle = true_anomaly * (2*np.pi / 360)
            r = (a * (1 - e**2)) / (1 + e * np.cos(satellite_angle))
            polar_x = r * np.cos(satellite_angle)
            polar_y = r * np.sin(satellite_angle)
            polar_z = (0) * satellite_angle
            pos_sat_orbit = np.matrix([polar_x, polar_y, polar_z]) # Current coordinate, no integration along orbit
        
        rot_matrix_mash = Minc * Masc * Mper

        # Newly oriented satellite trajectory
        oriented_sat_orbit = rot_matrix_mash * pos_sat_orbit.T
        oriented_sat_orbit = oriented_sat_orbit.T
        
        if satellite_body == False:
            # Plot Ellitical Orbit of Satellite 
            orb_x = sum(oriented_sat_orbit[:,0].tolist(), [])
            orb_y = sum(oriented_sat_orbit[:,1].tolist(), [])
            orb_z = sum(oriented_sat_orbit[:,2].tolist(), [])
        else: 
            sat = oriented_sat_orbit.flatten()
            orb_x = sat[0,0] # All correspond to position of physical satellite body at current epoch
            orb_y = sat[0,1] # extracted from space-track.org
            orb_z = sat[0,2]
        
        return pos_sat_orbit, orb_x, orb_y, orb_z
   
    def define_earth_grid(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1,projection = '3d')
        self.fig.set_size_inches(10,10)
        
        r_earth, x, y, z, polar_x_e, polar_y_e, polar_z_e = self.earth()
        self.ax.plot_wireframe(x, y, z, alpha = 0.2, color = 'blue', rstride = 4, cstride = 4) # Plot Cage of Earth
        self.ax.plot(polar_x_e, polar_y_e, polar_z_e, color = 'blue', linestyle = '-') # Plot Earth's Equator
        self.ax.set_axis_on()
        
    def plot(self, file, Name, index, show_name = False):
        
        theta = np.linspace(0, 2 * np.pi, 100)
        object_name, a, e, inc, asc, per, true_anomaly = self.extract_orbital_elements_sat(file, index)
        Minc, Masc, Mper = self.orientation_sat(inc, asc, per) 
        pos_sat_orbit, orb_x, orb_y, orb_z = self.orbit_pos(a, e, inc, asc, per, true_anomaly, satellite_body = False)
        pos_sat, satellitex, satellitey, satellitez = self.orbit_pos(a, e, inc, asc, per, true_anomaly, satellite_body = True)
        
        # LABEL AXIS for reference point
        # X ~ Vernal Equinox
        self.ax.plot([0,7500],[0,0],[0,0],'r:')
        self.ax.plot([8000],[0],[0],'r>')
        self.ax.text(8300,0,0,s='X (Vernal Equinox)', fontsize=12,color='black')

        # Create Y-axis Label
        self.ax.plot([0,0],[0,7500],[0,0],'r:')
        self.ax.plot([0],[8500],[0],'r>')
        self.ax.text(0,9400,0,s='Y',fontsize=12,color='black')

        # Create Z-axis Label
        self.ax.plot([0,0],[0,0],[0,7500],'r:')
        self.ax.plot([0],[0],[8000],'r^')
        self.ax.text(0,0,8300,s='Z', fontsize=12,color='black')
        
        # Set labels of grid
        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        self.ax.set_zlabel('Z (km)')
        
        # Plot Orbit Trajectory + Physical Body of Satellite
        self.ax.plot(orb_x, orb_y, orb_z, color = 'black', linestyle = '--', linewidth = 1.2) # Plot Satellite Orbit
        self.ax.plot([satellitex], [satellitey], [satellitez], 'go') # Plot body of satellite
        self.ax.plot([0, satellitex], [0, satellitey], [0, satellitez], 'g-') # Plot vector extending to satellite
        
        if show_name == True:
            self.ax.text(satellitex, satellitey, satellitez, s = object_name, fontsize = 9)
        else: 
            pass
        
        # Set tight_fit layout of 3d axis
        xyzlim = np.array([self.ax.get_xlim3d(), self.ax.get_ylim3d(), self.ax.get_zlim3d()]).T
        XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
        self.ax.set_xlim3d(XYZlim)
        self.ax.set_ylim3d(XYZlim)
        self.ax.set_zlim3d(XYZlim * 3/4)
        self.ax.set_title(f'Constellation Plot of {Name}', fontsize = 20)
        
                    
       
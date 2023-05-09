import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import numpy as np
from PIL import Image

def intro():
    import matplotlib.pyplot as plt
    import pandas as pd
    import streamlit as st
    import numpy as np
    from PIL import Image

    #st.set_page_config(page_title = 'ORBYTES')
    st.header('ORBYTES:')
    st.subheader('Orbit Representation BY Tracking Earth Satellites')
    st.sidebar.success('Select a demo above.')


    image = Image.open('Images/orbytes.jpeg')
    st.image(image, caption = 'DALL-E 2 generated image of ORBYTES')

    st.subheader('ORBYTES provides users the ability to interact with an API interface to :green[extract, visualize, and integrate] satellite trajectories actively tracked by space-track.org.')

    st.subheader('In addition to mapping the orbital paths of payloads in regions of LEO to GEO, ORBYTES functions as a platform for users to gain intuition about the beauty and complexity of orbital mechanics. ')

    st.subheader('This knowledge is established through step-by-step tutorials beginning with a lesson on orbital elements (differentiating between ones that :red[1) define the physical shape of an orbit] vs :red[2) how an orbit is oriented in space]) followed by a live demonstration of how current constellations of satellites assume their orbital trajectories. ')
                        
    st.subheader('Throughout this tutorial, the user will:' )

    st.write('* **Establish intuition about coordinate transformations within the Earth-Centered Inertial Frame (ECI)** ')
    st.write('* **Learn about the 7 main orbital elements that define the shape, structure, and orientation of an orbit**')
    st.write('* **Access archives of data on space-track.org using its Application Programming Interface (API)**')
    st.write('* **Modify the orbital architecture of orbital shells surrounding Earth to learn about the placement of satellites within larger mega-constellations**')
    st.write('* **Witness the evolution of the Mean Keplerian Elements of a satellite due to gravitational perturbations from those on adjacent orbital planes**')

def orbital_element_demo():
    import matplotlib.pyplot as plt
    import pandas as pd
    import streamlit as st
    import numpy as np
    from PIL import Image
    
    def eccentric_anomaly(file, index = 0, iterations = 500): # file is Pandas DataFrame
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        
        if mean_anomaly >= 0.6:
            init = np.pi
        else:
            init = mean_anomaly
        for i in range(iterations):
            ecc_anomaly = init - (init - eccentricity * np.sin(init) - mean_anomaly) / (1.0 - eccentricity * np.cos(init))
        
        return ecc_anomaly

    st.header('Orbital Elements Tutorial')

    st.write('This demo illustrates the basics of orbital mechanics. The user is able to experiment with the realationships that exist between the orbital parameters and maneuver how they change the trajectory of the satellite in question. We begin our introduction to 3D projections in this tutorial. ')

    image = Image.open('Images/orbit.png')
    st.image(image, caption = '2D Schematic of Orbital Element Relationships')

    st.write('Below you will find a list of orbital elements which help characterize the orbit of a satellite. Together with these parameters one is able to define the orbital plane on which a satellite will assume is respective trajectory around any central object (in this case, the Earth).')  
    st.write('**:red[Semi-Major Axis (a):]** Distance from the center of an ellipse to the most distant point on the ellipse.')
    st.write('**:red[Eccentricity (e):]** Ratio of the distances between center & focus and semi-major axis.')
    st.write('**:red[Inclination (i):]** The angle between the plane of reference (equatorial plane) and the orbital plane of the satellite. ')
    st.write('**:red[Argument of Perigee ($$\omega$$):]** The angle within the orbital plane of the satellite from the lines of nodes to the line of perigee.')
    st.write('**:red[Longitude of the Ascending Node ($$\Omega$$):]** The angle within the plane of reference (equatorial plane) between the vernal equinox and the line of nodes.') 
    st.write('**:red[True Anomaly (v):]** The angle within the orbital plane between the line of perigee and the location of satellite along its true trajectory.')

    st.write('**IMPORTANT POINT:** We need to determine how an orbit is oriented in space. Orientation is determined by 3 rotations. 1) About the Line of Notes (ref angle: Inclination), 2) About the Z-Axis (ref angle: longitude of the ascending node), and 3) Rotation along Angular Momentum vector (perp. orbital plane & ref angle: argument of perigee)')

    st.subheader('Select Orbital Parameters:')
    a = st.slider('**:blue[Semi-Major Axis (in km):]**', min_value = 6000, max_value = 26000, step = 1)
    e = st.slider('**:blue[Eccentricity :]**', min_value = 0.0, max_value = 1.0, step = 0.00001)
    i = st.slider('**:blue[Inclination (in degrees):]**', min_value = 0.0, max_value = 70.0, step = 0.5)
    per = st.slider('**:blue[Argument of Perigee (in degrees):]**', min_value = 0.0, max_value = 360.0, step = 0.5)
    asc = st.slider('**:blue[Longitude of the Ascending Node (in degrees):]**', min_value = 0.0, max_value = 360.0, step = 0.5)
    ano = st.slider('**:blue[True Anomaly (in degrees):]**', min_value = 0.5, max_value = 360.0, step = 0.5)


    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    fig.set_size_inches(10,10)

    # Plot Earth
    unit_vector = (1,1,1)
    r_x, r_y, r_z = [6378 / np.sqrt(component) for component in unit_vector]
    theta = np.linspace(0, 2*np.pi, 100)
    phi = np.linspace(0, np.pi, 100)
    x = r_x * np.outer(np.cos(theta), np.sin(phi))
    y = r_y * np.outer(np.sin(theta), np.sin(phi))
    z = r_z * np.outer(np.ones_like(theta), np.cos(phi)) 

    # Orientation components
    i = i / ((2*np.pi)/360) # in radians
    R = np.matrix([[1, 0, 0], # rotation about the x-axis (actually, line of nodes)
                    [0, np.cos(i), -np.sin(i)], 
                    [0, np.sin(i), np.cos(i)]])
    asc = asc / ((2*np.pi)/360) # in radians
    R2 = np.matrix([[np.cos(asc), -np.sin(asc), 0],
                    [np.sin(asc), np.cos(asc), 0],
                    [0, 0, 1]]) # rotation about the Z-axis 
    per = per / ((2*np.pi)/360) # in radians
    R3 = np.matrix([[np.cos(per), -np.sin(per), 0], 
                    [np.sin(per), np.cos(per), 0],
                    [0, 0, 1]]) # Rotation about the Z-axis (actually, angular momentum vector)

    # Polar equation 
    r = (a * (1 - e**2))/ (1 + e * np.cos(theta))
    polar_x = r * np.cos(theta)
    polar_y = r * np.sin(theta)
    polar_z = 0 * theta

    # Equatorial plane
    a_earth = 6378
    e_earth = 0
    r_earth = (a_earth * (1 - e_earth**2))/ (1 + e_earth * np.cos(theta))
    polar_x_earth = r_earth * np.cos(theta)
    polar_y_earth = r_earth * np.sin(theta)
    polar_z_earth = 0 * theta

    points_e = np.matrix(list(zip(polar_x_earth, polar_y_earth, polar_z_earth)))
    pts_e = points_e.T
    point_vector_e = pts_e.T
    x1_e, y1_e, z1_e = point_vector_e[:, 0].A.flatten(), point_vector_e[:, 1].A.flatten(), point_vector_e[:, 2].A.flatten()
    ax.plot(x1_e, y1_e, z1_e, color = 'blue', linestyle = '-')

    # Plot Satellite 
    sat = ano * ((2*np.pi)/360)
    satr = r = (a * (1 - e**2))/ (1 + e * np.cos(sat))
    satx = satr * np.cos(sat)
    saty = satr * np.sin(sat)
    satz = satr * 0

    satpoints = np.matrix([satx, saty, satz])
    sat_pts = (R * R2 * R3) * satpoints.T
    satellite_a = sat_pts.flatten()
    satellitex, satellitey, satellitez = satellite_a[0,0], satellite_a[0,1], satellite_a[0,2]
    ax.plot([0, satellitex], [0, satellitey], [0, satellitez], 'g-')
    ax.plot([satellitex],[satellitey], [satellitez], 'go')
    ax.text(satellitex, satellitey, satellitez, s = 'SATELLITE', fontsize = 11)

    # Plot Orbit 
    points = np.matrix(list(zip(polar_x, polar_y, polar_z)))
    pts = (R * R2 * R3) * points.T
    point_vector = pts.T
    x1, y1, z1 = point_vector[:, 0].A.flatten(), point_vector[:, 1].A.flatten(), point_vector[:, 2].A.flatten()
    ax.plot_wireframe(x,y,z, alpha = 0.2, color = 'blue', rstride = 4, cstride = 4)
    ax.plot(x1, y1, z1, color = 'black', linestyle = '-')

    # Create X-Axis Label
    ax.plot([0,7500],[0,0],[0,0],'r:')
    ax.plot([8000],[0],[0],'r>')
    ax.text(8300,0,0,s='X (Vernal Eqx)', fontsize=12 ,color='black')

    #Create Y-axis Label
    ax.plot([0,0],[0,7500],[0,0],'r:')
    ax.plot([0],[8500],[0],'r>')
    ax.text(0,9400,0,s='Y',fontsize=12,color='black')

    #Create Z-axis Label
    ax.plot([0,0],[0,0],[0,7500],'r:')
    ax.plot([0],[0],[8000],'r^')
    ax.text(0,0,8300,s='Z', fontsize=12,color='black')

    #Create Z-axis Label
    ax.plot([0,0],[0,0],[0,7500],'m-')
    ax.plot([0],[0],[7500],'m<')
    #ax.text(0,0,8000,s='axis', fontsize=15,color='b')

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),      
                       ax.get_zlim3d()]).T
    XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
    ax.set_xlim3d(XYZlim)
    ax.set_ylim3d(XYZlim)
    ax.set_zlim3d(XYZlim * 3/4)
    ax.set_title('Orbit Modeler for User Input Parameters')
    #ax.view_init(elev = 0, azim = 0)

    st.pyplot(fig)

def space_model():
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
    
    st.header('True Constellation Modeling using Space-Track.org Data')

    st.write('This demo illustrates the intermediate step in our efforts to comprehend orbital mechanics. Together with the knowledge of orbital trajectories garnered from the previous tutorial, the user is able to choose a constellation of satellites recently launched and plot their orbits around the central Earth. These satellites have a wide range of uses from telecommunications to scientific to military/government-operated missions. Once again we utilize 3D projections in an effort to simulate the most accurate orbital paths of payloads from LEO to GEO. ')

    st.write('**Satellites must possess this column data to be considered in orbit around Earth.**')
    st.write('* The value associated with the CENTER_NAME keyword shall be ‘EARTH’.')
    st.write('* The value associated with the REF_FRAME keyword shall be ‘TEME’ (see annex A).')
    st.write('* The value associated with the TIME_SYSTEM keyword shall be ‘UTC’.')

    image = Image.open('Images/orientation.png')
    st.image(image, caption = '2D Image of Orientation of Orbital Planes which satellites occupy')     

    # Need orbital elements extracted of specific satellites 
    option = st.selectbox('**:blue[Which Constellation of Satellites would you like to plot?]**',
                      ('STARLINK','STARLINK 1:20', 'GALILEO', 'GLONASS', 'GLONASS 1:20', 'GLOBALSTAR',
                      'GLOBALSTAR 1:20'))
    st.write('You selected:', option)

    def eccentric_anomaly(file, index = 0, iterations = 500): # file is Pandas DataFrame
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
            
    def select_data(output_file, **kwargs):
        reduced_data = pd.read_csv(output_file)
        for key, value in kwargs.items():
            reduced_data = reduced_data.loc[reduced_data['{}'.format(key)].str.contains('{}'.format(value), case = False)]
        return reduced_data.reset_index(drop = True)

    def extract_orbital_elements(file, index = 0):
        whole_file = file
        # Calls 'plot_orbit' function using the TLE elements defined in output file
        row_orbital_elements = file.loc[file.index[index]]
        object_name = file['OBJECT_NAME'].iloc[index]
        semi_major_axis = file['SEMIMAJOR_AXIS'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        inclination = file['INCLINATION'].iloc[index]
        right_ascension = file['RA_OF_ASC_NODE'].iloc[index]
        argument_pericenter = file['ARG_OF_PERICENTER'].iloc[index] # Angle between lines of nodes (ascending node) and periapsis 
        
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        ecc_anomaly = eccentric_anomaly(whole_file)
        true_anomaly = np.arctan((np.sqrt(1 - eccentricity**2) * np.sin(ecc_anomaly)) / (np.cos(ecc_anomaly) - eccentricity))# Solve for true anomaly using mean anomaly and eccentricity : defines the position of a body along the ellipse at a specific time ('epoch') - angle between periapsis & current location along orbit
    
        return object_name, semi_major_axis, eccentricity, inclination, right_ascension, argument_pericenter, true_anomaly

    full = pd.read_csv('Example_csv_set/all_data.csv')
    if option == 'STARLINK':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'STARLINK-*')
    elif option == 'STARLINK 1:20':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'STARLINK-*')
        data = data[0:20]
    elif option == 'GALILEO':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GALILEO-*')
    elif option == 'GLONASS':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLONASS-*')
    elif option == 'GLONASS 1:20':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLONASS-*')
        data = data[0:20]
    elif option == 'GLOBALSTAR':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLOBALSTAR-*')
    elif option == 'GLOBALSTAR 1:20':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLOBALSTAR-*')
        data = data[0:20]
    else: 
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = f'{option}-*')
        
    data.to_csv('new_data.csv')
    st.header(f'Sample Data for {option} constellation given by Space-Track.org')
    st.write(data.head())    

    # CREATE PLOTS
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection = '3d')
    fig.set_size_inches(10,10)

    # Plot Earth
    unit_vector = np.array([1,1,1])
    earth_pos_vector = []
    for cartesian_component in unit_vector:
        earth_pos_vector.append(6378 / cartesian_component)
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
    r_earth = (6378 * (1 - 0**2)) / (1 + 0 * np.cos(theta))
        
    # Polar coordinates of Earth's Equatorial Plane projection
    polar_x_e = r_earth * np.cos(theta)
    polar_y_e = r_earth * np.sin(theta)
    polar_z_e = r_earth * 0
        
    position_earth_3d = np.matrix(list(zip(polar_x_e, polar_y_e, polar_z_e)))

    # Plot Earth
    ax.plot_wireframe(x,y,z, alpha = 0.2, color = 'blue', rstride = 4, cstride = 4)
    ax.plot(polar_x_e, polar_y_e, polar_z_e, color = 'blue', linestyle = '-')
    ax.set_axis_on()


    for i in range(len(data)):
        object_name, a, e, inc, asc, per, true_anomaly = extract_orbital_elements('new_data.csv', index = i)
    
        # Orientation components
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
        # Polar equation for orbit
        r = (a * (1 - e**2))/ (1 + e * np.cos(theta))
        polar_x = r * np.cos(theta)
        polar_y = r * np.sin(theta)
        polar_z = 0 * theta
        pos_sat_orbit = np.matrix(list(zip(polar_x, polar_y, polar_z)))
        
        rot_matrix_mash = Minc * Masc * Mper

        # Newly oriented satellite trajectory
        oriented_sat_orbit = rot_matrix_mash * pos_sat_orbit.T
        oriented_sat_orbit = oriented_sat_orbit.T
        
        orb_x = sum(oriented_sat_orbit[:,0].tolist(), [])
        orb_y = sum(oriented_sat_orbit[:,1].tolist(), [])
        orb_z = sum(oriented_sat_orbit[:,2].tolist(), [])
        
        # Plot Satellite Orbit
        ax.plot(orb_x, orb_y, orb_z, color = 'black', linestyle = '--', linewidth = 1.2) # Plot Satellite Orbit
        
        # Plot body of satellite (not orbit)
        satellite_angle = true_anomaly * (2*np.pi / 360)
        r = (a * (1 - e**2)) / (1 + e * np.cos(satellite_angle))
        polar_x = r * np.cos(satellite_angle)
        polar_y = r * np.sin(satellite_angle)
        polar_z = (0) * satellite_angle
        pos_sat_orbit = np.matrix([polar_x, polar_y, polar_z]) # Current coordinate, no integration along orbit
        
        rot_matrix_mash = Minc * Masc * Mper
        oriented_sat_orbit = rot_matrix_mash * pos_sat_orbit.T
        oriented_sat_orbit = oriented_sat_orbit.T
        
        sat = oriented_sat_orbit.flatten()
        satellitex = sat[0,0] # All correspond to position of physical satellite body at current epoch
        satellitey = sat[0,1] # extracted from space-track.org
        satellitez = sat[0,2]
        
        # Plot body of satellite
        ax.plot([0, satellitex], [0, satellitey], [0, satellitez], 'g-')
        ax.plot([satellitex],[satellitey], [satellitez], 'go')
        ax.text(satellitex, satellitey, satellitez, s = object_name, fontsize = 9)

    # Create X-Axis Label
    ax.plot([0,7500],[0,0],[0,0],'r:')
    ax.plot([8000],[0],[0],'r>')
    ax.text(8300,0,0,s='X (Vernal Eqx)', fontsize=12 ,color='black')

    #Create Y-axis Label
    ax.plot([0,0],[0,7500],[0,0],'r:')
    ax.plot([0],[8500],[0],'r>')
    ax.text(0,9400,0,s='Y',fontsize=12,color='black')

    #Create Z-axis Label
    ax.plot([0,0],[0,0],[0,7500],'r:')
    ax.plot([0],[0],[8000],'r^')
    ax.text(0,0,8300,s='Z', fontsize=12,color='black')

    # Set Grid Labels
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # Set tight_fit layout of 3d axis
    xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),      
                       ax.get_zlim3d()]).T
    XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
    ax.set_xlim3d(XYZlim)
    ax.set_ylim3d(XYZlim)
    ax.set_zlim3d(XYZlim * 3/4)
    ax.set_title(f'Constellation Plot of {option} Satellites', fontsize = 20)
    #ax.view_init(elev = 0, azim = 0)

    st.pyplot(fig)

def rebound_demo():
    import matplotlib.pyplot as plt
    import pandas as pd
    import streamlit as st
    import numpy as np
    from PIL import Image
    import numpy as np
    import rebound
    from mpl_toolkits.mplot3d import Axes3D
    import mpld3
    import streamlit.components.v1 as components
    import plotly.express as px
    
    
    st.header('Evolution of Satellite Orbital Elements (specifically a, e, and inc)')

    st.write('In this demo we begin to analyze the evolution of each satellites orbital elements as a direct consequence of the proximity as well as gravitational perturbations that ensue. We once again refine our approach and understanding of orbital mechanics. We are also able to visualize through the use of high-order symplectic integrators the evolution in three important orbital elements amongst each satellite. ')
    
    def select_data(output_file, **kwargs):
        reduced_data = pd.read_csv(output_file)
        for key, value in kwargs.items():
            reduced_data = reduced_data.loc[reduced_data['{}'.format(key)].str.contains('{}'.format(value), case = False)]
        return reduced_data.reset_index(drop = True)
    
    def eccentric_anomaly(file, index = 0, iterations = 500): # file is Pandas DataFrame
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        
        if mean_anomaly >= 0.6:
            init = np.pi
        else:
            init = mean_anomaly
        for i in range(iterations):
            ecc_anomaly = init - (init - eccentricity * np.sin(init) - mean_anomaly) / (1.0 - eccentricity * np.cos(init))
        
        return ecc_anomaly
    
    def plot(a,b,c,d, object_name):
        import matplotlib.pyplot as pyplot

        rows, cols = 3, 1
        fig, (ax1, ax2, ax3) = pyplot.subplots(rows, cols, sharex = True)
        ax1.plot(a/(2*np.pi), c)
        ax1.set_ylabel('Semi-major Axis', fontsize = 8)
        ax1.tick_params(labelsize = 8)
        ax1.set_xlim(0,1000)
        
        ax2.plot(a/(2*np.pi), b)
        ax2.set_ylabel('Eccentricity', fontsize = 8)
        ax2.tick_params(labelsize = 7.5)
        ax2.set_xlim(0,1000)
        
        ax3.plot(a/(2*np.pi), d * (360/(2*np.pi)))
        ax3.set_ylabel('Inclination', fontsize = 8)
        ax3.tick_params(labelsize = 8)
        ax3.set_xlabel('Time (years)')
        ax3.set_xlim(0, 1000)

        ax1.set_title(f'Orbital Evolution of Satellite {object_name}', fontsize = 10, pad = 20)
        st.pyplot(fig)

    def extract_orbital_elements(file, index = 0):
        # Calls 'plot_orbit' function using the TLE elements defined in output file
        row_orbital_elements = file.loc[file.index[index]]
        object_name = file['OBJECT_NAME'].iloc[index]
        semi_major_axis = file['SEMIMAJOR_AXIS'].iloc[index]
        eccentricity = file['ECCENTRICITY'].iloc[index]
        inclination = file['INCLINATION'].iloc[index]
        right_ascension = file['RA_OF_ASC_NODE'].iloc[index]
        argument_pericenter = file['ARG_OF_PERICENTER'].iloc[index] # Angle between lines of nodes (ascending node) and periapsis 
        
        mean_anomaly = file['MEAN_ANOMALY'].iloc[index]
        ecc_anomaly = eccentric_anomaly(file)
        true_anomaly = np.arctan((np.sqrt(1 - eccentricity**2) * np.sin(ecc_anomaly)) / (np.cos(ecc_anomaly) - eccentricity))# Solve for true anomaly using mean anomaly and eccentricity : defines the position of a body along the ellipse at a specific time ('epoch') - angle between periapsis & current location along orbit
    
        return object_name, semi_major_axis, eccentricity, inclination, right_ascension, argument_pericenter, mean_anomaly

    # Need orbital elements extracted of specific satellites 
    option = st.selectbox('**:blue[Which Constellation of Satellites would you like to plot?]**',
                          ('STARLINK','STARLINK 1:20', 'GALILEO', 'GLONASS', 'GLONASS 1:20'))
    st.write('You selected:', option)

    full = pd.read_csv('Example_csv_set/all_data.csv')
    if option == 'STARLINK':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'STARLINK-*')
    elif option == 'STARLINK 1:20':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'STARLINK-*')
        data = data[0:20]
    elif option == 'GALILEO':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GALILEO-*')
    elif option == 'GLONASS':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLONASS-*')
    elif option == 'GLONASS 1:20':
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = 'GLONASS-*')
        data = data[0:20]
    else: 
        data = select_data('Example_csv_set/all_data.csv', OBJECT_NAME = f'{option}-*')
    
    index = st.slider('**:blue[SATELLITE INDEX]**', min_value = 0, max_value = len(data), step = 1)

    sim = rebound.Simulation()
    sim.add('Earth')
    for i in range(len(data)):
        object_name, a, e, i, asc, per, mean = extract_orbital_elements(data, index = i)
        f = rebound.M_to_f(e, mean)
        sim.add(a = a, e = e, inc = np.radians(i), Omega = np.radians(asc), omega = np.radians(per), f = f)
    
    particles = sim.particles
    sim.integrator = 'janus'

    torb = 2. *np.pi
    times = np.linspace(0*torb, 1000*torb, 200)
    sim.dt = 1e-3

    a = np.zeros(200)
    e = np.zeros(200)
    inc = np.zeros(200)
    
    for i, time in enumerate(times):
        sim.integrate(time, exact_finish_time=0)
        orbit = sim.particles[index].calculate_orbit(primary = sim.particles[0])
        a[i] = orbit.a
        e[i] = orbit.e
        inc[i] = orbit.inc

    orbital_prop_SS = np.vstack((times, e, a, inc))
    orbital_evol = np.transpose(orbital_prop_SS)
    df = pd.DataFrame(orbital_evol, columns = ['times(yr)', 'e', 'a', 'inc'])
    object_name_final, a, e, i, asc, per, mean = extract_orbital_elements(data, index = index)
    
    plot(df['times(yr)'], 1- df['e'], df['a'], df['inc'], object_name_final)
    

page_names_to_funcs = {
    "—": intro,
    "Orbital Elements Tutorial": orbital_element_demo,
    "Space-Track.org Data Visualization": space_model,
    "Evolution of Satellite Orbital Element": rebound_demo
}

demo_name = st.sidebar.selectbox("Choose a demo", page_names_to_funcs.keys())
page_names_to_funcs[demo_name]()
    
    
    
    
    
    
    
    
    
    
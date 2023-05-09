import os
import json
import configparser
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class setup_files():
    """ Class that has 5 methods for initialization of data by allowing user to access data
    on space-track.org using its Application Programming Interface (API)
    
    Attributes
    __________
    name_ini: str
        File name chosen for INI File (where login credentials (i.e. username and password) for access to space-track.org are held)
    username: str
        Personal Username for space-track.org
    password: str
        Personal Password for space-track.org
    output_file_name: str
        Name chosen for excel (or .csv) file which contains extracted data in GP format
    
    
    Methods
    _______
    create_ini()
        Create INI file, a configuration file that incorporates key-value pairs to initalize session,
        which contains user's Username & Password to access space-track.org's
        database, as well as the name of the output file to write extracted data to
    access_ini()
        Access INI file to obtain necessary credentials for space-track.org
    generate_query()
        Implement very own predicates when generating query for extacting data from space-track.org 
        (different columns for filtering predicates can be found at 
        https://www.space-track.org/basicspacedata/modeldef/class/gp/format/html) 
    url_request()
        User input commands to access satellite data via url
    satellite_csv()
        Create '.csv' file of each satellite's 1) Mean Keplerian Elements & 2) TLE Related Parameters
    extract_csv()
        Run methods in order, ultimately producing a '.csv' file comprised of satellite data (chosen by user)
    
    """
    def __init__(self, name_ini, username, password):
        # Initialize INI file components
        self.name_ini = name_ini +'.ini'
        self.username = username
        self.password = password
        
    def create_ini(self):
        # Create INI file and save to working directory
        new_ini = open(self.name_ini, 'w')
        new_ini.write('[Space-Track credentials]\n')
        new_ini.write(f'Username = {self.username}\n')
        new_ini.write(f'Password = {self.password}\n')
        new_ini.close()
        
    def access_ini(self):
        # Extract contents of INI file which will be used to log onto space-track.org
        config = configparser.ConfigParser()
        config.read(f'./{self.name_ini}')
        user = config['Space-Track credentials']['Username']
        pasw = config['Space-Track credentials']['Password']
        sign_in = {'identity': user, 'password': pasw}
        return sign_in
    
    def generate_query(self, generate = False, orderby, **predicates):
        # MUST ALWAYS MENTION AN 'ORDERBY' ARGUMENT (i.e. OBJECT_NAME, OBJECT_ID, etc.)
        # If generating own query to extract specific category of data, enter all predicates/constraints on constellation to form url 
        query_str = ""
        if generate == True:
            # Use general perturbations class (gp) - efficient collection of keplerian orbital elements
            # for each earth-orbiting object tracked by 18th Space Defense Squadron
            class_type = 'class/gp/' 
            query_str += class_type
    
            for key, value in predicates.items():
                if key == 'orderby':
                    orderby_val = value
                    pass
                else:
                    query_str += f'{key}/{value}/'
                    
            query_str += f'orderby/{orderby_val}/format/json'
        else:
            pass
        
        return query_str
    
    def url_request(self, all_data = False, payload = False, manmade = False, **kwargs):
        urlBase = 'https://www.space-track.org/'
        RequestLogin = 'ajaxauth/login/'
        
        # 'basicspacedata': data available to all members with valid user & pass 
        # 'query': action used to extract data specified with the rest of the URL parameters
        RequestControllerAction = 'basicspacedata/query/'  
        
        if all_data == True: # Extract all data available on space-track.org (Payloads, debris, rocket parts, etc.) (Download time ~ 30 min)
            # Retrieve newest propagable element set for all objects which are currently on their assigned trajectories
            REQDATA = 'class/gp/decay_date/null-val/epoch/%3Enow-30/orderby/norad_cat_id/format/json'
        
        elif payload == True:  # Extract all satellite data that has not yet decayed (Download time ~ 10 min)
            REQDATA = 'class/gp/decay_date/null-val/epoch/%3Enow-30/REF_FRAME/TEME/CENTER_NAME/EARTH/TIME_SYSTEM/UTC/OBJECT_TYPE/PAYLOAD/orderby/norad_cat_id/format/json'
        
        elif manmade == True:
            query_str = self.generate_query(generate = True, **kwargs)
            REQDATA = query_str # Query generated from previous method
        
        return REQDATA
    
    def satellite_csv(self, sign_in, REQDATA, read_out_file):
        urlBase = 'https://www.space-track.org/'
        RequestLogin = 'ajaxauth/login/'
        RequestControllerAction = 'basicspacedata/query/'
        RequestData = REQDATA
        
        # use requests package to drive the RESTful session with space-track.org
        with requests.Session() as session: # session is object made with requests
        # run session in a with block to force session to close if we exit
    
            # 1) Log-in to space-track.org 
            req = session.post(urlBase + RequestLogin, data = sign_in)
    
            # 2) Input in-query for data extraction
            req = session.get(urlBase + RequestControllerAction + RequestData)
    
            # 3) Use json package to convert file format from json ==> list of dictionaries for Python to easily interpret
            DATA = json.loads(req.text) # JavaScript Object Notation (JSON) format used to send data as text == objects are defined by their name & value -- similar to python dictionaries

    
            for i in range(len(DATA)):
                if os.path.isfile(f'{read_out_file}.csv'):
                    df_read = pd.read_csv(f'{read_out_file}.csv')
                else:
                    df = pd.DataFrame(list(), columns = ("CREATION_DATE", "OBJECT_NAME", "OBJECT_ID", "CENTER_NAME",
                                                          "REF_FRAME", "TIME_SYSTEM", 'MEAN_ELEMENT_THEORY', 'EPOCH', 
                                                          'MEAN_MOTION', 'ECCENTRICITY', 'INCLINATION', 'RA_OF_ASC_NODE',
                                                          'ARG_OF_PERICENTER', 'MEAN_ANOMALY', 'EPHEMERIS_TYPE', 'REV_AT_EPOCH',
                                                          'BSTAR', 'MEAN_MOTION_DOT', 'MEAN_MOTION_DDOT', 'SEMIMAJOR_AXIS', 
                                                          'PERIOD', 'APOAPSIS', 'PERIAPSIS','OBJECT_TYPE', 
                                                          'LAUNCH_DATE', 'DECAY_DATE'))
                    df.to_csv(f'{read_out_file}.csv')
                    df_read = pd.read_csv(f'{read_out_file}.csv')
               
                # Turn into Function -> Table
                creation_date = DATA[i]['CREATION_DATE']
                obj_name = DATA[i]['OBJECT_NAME']
                obj_id = DATA[i]['OBJECT_ID']
                center_name = DATA[i]['CENTER_NAME']
                ref_name = DATA[i]['REF_FRAME'] 
                time_name = DATA[i]["TIME_SYSTEM"]
                mean_el_th = DATA[i]["MEAN_ELEMENT_THEORY"]
                epoch = DATA[i]['EPOCH']
                mean_motion = DATA[i]['MEAN_MOTION']
                ecc = DATA[i]['ECCENTRICITY'] 
                inc = DATA[i]["INCLINATION"]
                asc_node = DATA[i]["RA_OF_ASC_NODE"] 
                arg_pericenter = DATA[i]['ARG_OF_PERICENTER']
                mean_anomaly = DATA[i]['MEAN_ANOMALY']
                eph_type = DATA[i]['EPHEMERIS_TYPE']
                rev_epoch = DATA[i]['REV_AT_EPOCH']
                drag = DATA[i]['BSTAR']
                mean_motion_dot = DATA[i]['MEAN_MOTION_DOT']
                mean_motion_ddot = DATA[i]['MEAN_MOTION_DDOT'] 
                semi_axis = DATA[i]["SEMIMAJOR_AXIS"]
                period = DATA[i]["PERIOD"]
                apoapsis = DATA[i]['APOAPSIS']
                periapsis = DATA[i]['PERIAPSIS']
                obj_type = DATA[i]['OBJECT_TYPE'] 
                launch_date = DATA[i]["LAUNCH_DATE"]
                decay_date = DATA[i]["DECAY_DATE"] 

                if df_read.empty:
                    data = []
                    data.append([creation_date, obj_name, obj_id, center_name, ref_name, time_name, mean_el_th,
                                 epoch, mean_motion, ecc, inc, asc_node, arg_pericenter, mean_anomaly, eph_type, 
                                 rev_epoch, drag, mean_motion_dot, mean_motion_ddot, semi_axis, period, apoapsis, 
                                 periapsis, obj_type, launch_date, decay_date])
                    df = pd.DataFrame(data, columns = ("CREATION_DATE", "OBJECT_NAME", "OBJECT_ID", "CENTER_NAME",
                                                      "REF_FRAME", "TIME_SYSTEM", 'MEAN_ELEMENT_THEORY', 'EPOCH', 
                                                      'MEAN_MOTION', 'ECCENTRICITY', 'INCLINATION', 'RA_OF_ASC_NODE',
                                                      'ARG_OF_PERICENTER', 'MEAN_ANOMALY', 'EPHEMERIS_TYPE', 'REV_AT_EPOCH',
                                                      'BSTAR', 'MEAN_MOTION_DOT', 'MEAN_MOTION_DDOT', 'SEMIMAJOR_AXIS', 
                                                      'PERIOD', 'APOAPSIS', 'PERIAPSIS','OBJECT_TYPE', 
                                                      'LAUNCH_DATE', 'DECAY_DATE'))
                    df.to_csv(f'{read_out_file}.csv', index = False)

                else:
                    data = []
                    data.append([creation_date, obj_name, obj_id, center_name, ref_name, time_name, mean_el_th,
                                    epoch, mean_motion, ecc, inc, asc_node, arg_pericenter, mean_anomaly, eph_type, 
                                    rev_epoch, drag, mean_motion_dot, mean_motion_ddot, semi_axis, period, apoapsis, 
                                    periapsis, obj_type, launch_date, decay_date])
                    df = pd.DataFrame(data, columns = ("CREATION_DATE", "OBJECT_NAME", "OBJECT_ID", "CENTER_NAME",
                                                        "REF_FRAME", "TIME_SYSTEM", 'MEAN_ELEMENT_THEORY', 'EPOCH', 
                                                        'MEAN_MOTION', 'ECCENTRICITY', 'INCLINATION', 'RA_OF_ASC_NODE',
                                                        'ARG_OF_PERICENTER', 'MEAN_ANOMALY', 'EPHEMERIS_TYPE', 'REV_AT_EPOCH',
                                                        'BSTAR', 'MEAN_MOTION_DOT', 'MEAN_MOTION_DDOT', 'SEMIMAJOR_AXIS', 
                                                        'PERIOD', 'APOAPSIS', 'PERIAPSIS','OBJECT_TYPE', 
                                                        'LAUNCH_DATE', 'DECAY_DATE'))
                    df_read = pd.read_csv(f'{read_out_file}.csv')
                    df_read = pd.concat([df_read, df], ignore_index = True, axis = 0)
                    df_read.to_csv(f'{read_out_file}.csv', index = False)
                    
    def extract_csv(self, read_out_file, all_data = False, payload = False, generate = False, **kwargs):
        self.create_ini()
        sign_in = self.access_ini()

        if generate == True:
            DATA = self.url_request(manmade = True, **kwargs)
        elif payload == True:
            DATA = self.url_request(all_data = False, payload = True, manmade = False)
        else:
            DATA = self.url_request(all_data = True, payload = False, manmade = False)
        
        self.satellite_csv(sign_in = sign_in, REQDATA = DATA, read_out_file = read_out_file)
        

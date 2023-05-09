import os
import json
import configparser
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class reduce_data():
    """ Class that has 1 method for modification of the dataset from which the user will be crafting Orbit Models 
    Goal: If user has downloaded full, up-to-date file of 'all_data' or 'payloads' - gives them ability to create smaller DataFrame of objects pertinent to their current work
    
    Attributes
    __________
    output_file: str
        Name of the '.csv' file which contains all extracted satellite data
    
    Methods
    _______
    range_str()
        Modify table of data to only contain satellites of interest using string key word arguments (i.e. OBJECT_NAME, OBJECT_TYPE, EPOCH...)
    one_sided_range()
        Modify table of data to only contain satellites which are contained with respect to boundary condition (ex. Starlink satellite w/ ECCENTRICITY 'less than or equal to' 0.00058)
    two_sided_range()
        Modify table of data to only contain satellites which are contained in RANGE of specific key word argument (ex. Starlink Satellite w/ SEMIMAJOR_AXIS between 6500km and 6675km)
        
    
    """
    def __init__(self, output_file):
        self.output_file = output_file
        
    def range_str(self, **kwargs):
        # Resize data using str loc functions
        reduced_data = pd.read_csv(self.output_file)
        for key, value in kwargs.items():
            reduced_data = reduced_data.loc[reduced_data['{}'.format(key)].str.contains('{}'.format(value), case = False)]
        return reduced_data
    
    def one_sided_range(self, sign, **kwargs):
        # Resize data with one-sided range 
        reduced_data = pd.read_csv(self.output_file)
        for key, value in kwargs.items():
            if sign == 'greater than':
                reduced_data = reduced_data.loc[reduced_data[f'{key}'] > float(f'{value}')]
            if sign == 'greater than or equal to':
                reduced_data = reduced_data.loc[reduced_data[f'{key}'] >= float(f'{value}')]
            if sign == 'less than':
                reduced_data = reduced_data.loc[reduced_data[f'{key}'] < float(f'{value}')]
            if sign == 'less than or equal to':
                reduced_data = reduced_data.loc[reduced_data[f'{key}'] <= float(f'{value}')]
        return reduced_data
    
    def two_sided_range(self, greater_than, less_than, **kwargs): # Data that is DESIRED satisfies two conditions
        reduced_data = pd.read_csv(self.output_file)
        for key, value in kwargs.items():
            index_true = reduced_data[(reduced_data[f'{value}'] <= greater_than) | (reduced_data[f'{value}'] >= less_than)].index
            reduced_data.drop(index_true, inplace = True)
        return reduced_data
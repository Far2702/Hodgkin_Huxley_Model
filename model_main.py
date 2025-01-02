import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

class Hodgkin_Huxley_Model:
    
    def __init__(self,Resting_Potential=-65,External_Current=10,time_span=100,time_steps=30000):     # There is some issues for time_step less than 9000  
        """  
        The simulation starts from the normal resting membrane potential, but plotting begins from the user-defined potential(if given) to maintain constant reversal potentials despite changes in resting membrane potential.
        """
        self.Resting_potenital=Resting_Potential
        self.External_Current=External_Current
        self.time_span=time_span
        self.time_steps=time_steps
    
    def display_info(self):
        print("Resting Potential: ",self.Resting_potenital ,"mV")
        print("External Current: ",self.External_Current, 'mA')
        print("Time Span: ",self.time_span, 'ms')
        print("Time Steps: ",self.time_steps)
        
    def default_values(self):
        self.membrane_capacitance= 1 # μF/cm^2
        self.g_Na= 120 # mS/cm^2
        self.g_K= 36  # ms/cm^2
        self.g_l=0.3  # ms/cm^2
        self.E_Na= 50 # mV
        self.E_K = -77.0 # mV
        self.E_L = -54.5 # mV
    

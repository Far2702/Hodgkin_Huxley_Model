import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

class Hodgkin_Huxley_Model:
    
    def __init__(self,Resting_Potential=-65,External_Current=10,time_span=100,time_steps=10000):     # There is some issues for time_step less than 7000   
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
    
    def Sodium_conductance(self,t,V,h):     # Simplified Assumption I
        alpha_m=0.1*(V+40)/(1-np.exp(-(V+40)/10))   # m is the probability that the activation gate of Na+ is open
        beta_m=4*np.e**(-(V+65)/18)
        alpha_h=0.07*np.e**(-(V+65)/20)             # h is the probability that the inactivation gate of Na+ is open
        beta_h=1/(1+np.e**(-(V+35)/10))
        
        
        # if t==0:
        #     return 0
        # elif t>10 and t<=11.5:
        #     return self.g_Na
        # elif t>=11.5:
        #     return 0
        # else:
        #     return 0
        
    def Potassium_conductance(self,t,V,h):  
        
        if V==-55:
            alpha_n=0.1
        else:
            alpha_n=(0.01*(V+55))/(1-np.exp(-(V+55)/10))     # n is the probability that an individual subunit of K+ is open 
        beta_n= 0.125*np.e**(-(V+65)/80)
        n=alpha_n/(alpha_n+beta_n)
        n_new=n
        t_new=t
        t_new=t_new+h
        n_new=n_new+ h*((alpha_n*(1-n_new))-(beta_n*n_new))
        print(t_new)
        print(n_new)
        return self.g_K*((n_new)**4),t_new
        # if t==0:
        #     return 0
        # elif t>=12 and t<=40:
        #     return self.g_K
        # elif t>=40:
        #     return 0
        # else:
        #     return 0
          
    def plot_Membrane_Pot_vs_time(self):       # Solving differential equation using Euler Method
        V_m=[]
        t_0=0
        t=[]
        Vm_0=self.Resting_potenital
        V_new=Vm_0
        t_new=t_0
        h= self.time_span / self.time_steps
        
        while (t_new<=self.time_span):   # t_0=0 and Vm_0=-65mV
            V_m.append(V_new)
            t.append(t_new)
            K,t_new=self.Potassium_conductance(t_new,V_new,h)
            V_new=V_new+h*(-((self.Sodium_conductance(t_new)*(V_new-self.E_Na))/self.membrane_capacitance)-((K*(V_new-self.E_K))/self.membrane_capacitance)-((self.g_l*(V_new-self.E_L))/self.membrane_capacitance)+self.External_Current)
        
        fig,ax=plt.subplots()
        ax.set_xlim([min(t),max(t)])
        ax.set_ylim([min(V_m)-10,max(V_m)+10])
        
        animated_plot, = ax.plot([],[])        
        def update(frame):
            idx = frame*5    # Skipping 5 frames for fast animation
            animated_plot.set_data(t[:idx],V_m[:idx])
            return animated_plot,
        
        animation = FuncAnimation(fig=fig,func=update,frames=self.time_steps,interval=0,blit=True)

        plt.xlabel("Time(ms)")
        plt.ylabel("Potential(mV)")
        plt.show()
        print(V_m[:5])
        print(len(t))

##############--------------------------------------------
model=Hodgkin_Huxley_Model(External_Current=0)
model.default_values()
model.display_info()
model.plot_Membrane_Pot_vs_time()  
                
            
            
            
              
        

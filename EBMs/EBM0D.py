import numpy as np
import matplotlib.pyplot as plt






#Tune tau
S0=1361. # W/m^2 Solar constant 
alpha=0.3  # Albedo
T_goal=288. # K Temperature for tuning
sigma=5.67e-8 # W/m^2/K^4
# Tune tau
tau= S0/4 * (1-alpha) / (sigma * T_goal**4)



#Boltz
T0_B=(S0/4 * (1-alpha) / (sigma*tau))**0.25
P_hum_glob=0.034 # W/m^2 AHF today
dT_P_B=((S0/4*(1-alpha)+ P_hum_glob)/(tau*sigma))**0.25- T0_B






#LINEAR
# for (Mcguffie and Hendersson-Sellers, 2014)
A_M=-388.74 
B_M=2.17
T0_M=(S0/4*(1-alpha)-A_M)/B_M
dT_P_M=(S0/4*(1-alpha)-A_M+P_hum_glob)/B_M - T0_M
# for (Caldeira and Kasting, 1992)
def xi(co2):	# needs CO2 concentration in ppm
	return np.log(co2*1./300.)
def A_cald(co2):
	c=xi(co2)
	return -326.4 + 9.161 *c - 3.164 *c**2 + 0.5468* c**3
def B_cald(co2):
	c=xi(co2)
	return 1.953 - 0.04866 *c +0.01309* c**2 - 0.002577* c**3
CO2=388
A_C=A_cald(CO2) 
B_C=B_cald(CO2)
T0_C=(S0/4*(1-alpha)-A_C)/B_C
dT_P_C=(S0/4*(1-alpha)-A_C+P_hum_glob)/B_C - T0_C
print("Boltzmann: ", dT_P_B, "Linear Mcguffie and Henderson-Sellers (2014):", dT_P_M,  "Linear Caldeira and Kastings(1992):", dT_P_C)



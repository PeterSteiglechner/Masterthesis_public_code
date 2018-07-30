# Peter Steiglechner
# 30.7.2018
# This code reproduces the result for the EBM-0D with the OLR parametrisation of the conceptual model in Wilson and Gea-Banacloche (2012)

import numpy as np
import matplotlib.pyplot as plt
import sys

# constants and Parameter
v_0=667.5			#1/cm	# wavenumber with best absorption
cross_0=3.71e-23	#m^2	#crosssection peak of CO2 spectrum at v0
r_plus=0.086		#cm		#decrease in spectrum for v>v0
r_minus=0.092		#cm		#decrease in spectrum for v<v0
		
#for calculation of concentration of co2 molecules in atmosphere
N_A=6.022140857e23	# 1/mol	#Avogadro
p_atm=1013.		#mbar		#pressure of atm
R=8.31446	#J/mol/K		#Gas constant

#thickness of the atmosphere: see Wilson (is this a tuning parameter?)
L=8e3			#m

#Energy constants	#eV to make numbers larger/smaller
h=4.14e-15		#eVs	#planck const
c=3e8		#m/s		#speed of light
k=8.62e-5	#eV/K		#Boltzmann const

T_T=217.	#K	#Temp of Tropopause
		
xi_1=1-np.exp(-11000./L)	#Xi (measure for height) at Tropopause (11km)
xi_avg=(1-xi_1/2)			# average xi
    
v_bottom=10		#start of v array
v_top=4000		#end of v array
No=3000		# number of calculations
v_arr=np.linspace(v_bottom, v_top, num=No)
         
co2_today=380.
T_S_today=288.
	
alpha=0.3
S0=1361.
sigma=5.67e-8

##############################################################
#########     HELP FUNCTIONS             #####################
##############################################################

#returns number of CO2 molecules per m^3 near surface for a given ppm-Concentration and Temperature near the surface
def n0_calc( co2, T_S):		#in ppm and K
  return N_A*0.1*p_atm*co2/R/T_S*1e-3

#calculates cross section for an array of wavenumbers, for single wavenumber [1/cm] modify function (see Week4)
def cross( v):
  cross_section=np.empty([len(v)])
  for i in range(0,len(v)):
    if v[i]>v_0:
      r=r_plus
    else: 
      r=r_minus
    cross_section[i]=cross_0 * np.exp(-r* abs(v[i]-v_0))
  return cross_section

#calculates number of random walk steps to escape the atmosphere
def N( v, co2, T_S):	#v in 1/cm , co2 in ppm, T_S in K
  return (cross(v)*n0_calc(co2, T_S)*L)

#calulates the Boltzmannspectrum for a wavenumber[1/cm] array and T of emitter
def B( nu, T):	#nu in 1a/cm and T of emitter in K
  v_in_pro_m=nu*100	#in 1/m	
  v_in_pro_s=c*v_in_pro_m	#in 1/s
  spectrum=2*np.pi*h *v_in_pro_s**3 /c**2 *planck(v_in_pro_s, T)	#in eVs*1/s^3 / (m^2 /s^2 )= eV/m^2
  return spectrum*1.60218e-19	*v_in_pro_s/nu		#J/m^2 *1/s *1/(1/cm) = W/m^2 /cm^-1

#calculates the planck part of the Boltzmann spectrum given v in [1/s]
def planck( v_s, T): # v in [1/s]
  return 1./(np.exp(h*v_s/(k*T))-1.)

#multiplies Planck spectrum with Absorption due to CO2
def B_out( v, T_S, T_T, xi_avg, co2):	# v in 1/cm , co2 in ppm
  #T_S is emitting Temp, T_T is temperature of Tropopause at 11km (xi_1)
  # xi_avg is a measure for average height
  return B(v,T_S)*np.exp(-N(v, co2, T_S)*xi_avg) + (1-np.exp(-N(v, co2, T_S)*xi_avg) )*B(v, T_T) 

def integrate( B, v_arr):
  dv=v_arr[1:]-v_arr[:-1]
  return sum(abs(B[:-1] + B[1:])*0.5*dv)

def calc_x_non_co2( x_co2, T_S):
  x_tot=- 255.**4/ T_S**4  +1
  x_non_CO2=x_tot-x_co2
  return x_non_CO2

def calc_x_co2( co2, T_S):
  B_with_co2=B_out(v_arr, T_S, T_T, xi_avg, co2)	#today
  B_no_atm=B(v_arr,T_S)
  x_co2=integrate(B_with_co2 - B_no_atm, v_arr) / integrate(B_no_atm,v_arr)
  return x_co2
  
  

##############################################################
#########     IMPORTANT FUNCTIONS        #####################
##############################################################

def calc_OLR( T_S=T_S_today, co2=380 , co2_today=co2_today):
  #### x_non_co2 is calculated via the fix values for T_S, CO2, ...
  x_co2=calc_x_co2(co2_today, T_S)
  x_non_co2=calc_x_non_co2(x_co2, T_S)
  #print("T_S", T_S)
  #print("co2", co2)
  B_no_atm=B(v_arr,T_S)
  B_with_co2=B_out(v_arr, T_S, T_T, xi_avg, co2)
  return integrate((B_with_co2- x_non_co2 * B_no_atm), v_arr)


# for single P_hum and co2 values. 	
def calc_T(co2=co2_today, P_hum=0.034):
  T=[T_S_today]
  h=0.01
  tolerance=0.0001
  x_co2_today=calc_x_co2(co2_today, T_S_today)
  x_non_co2=calc_x_non_co2(x_co2_today, T_S_today)
  while True:
    x_co2=calc_x_co2(co2, T[-1])
    T_new=T[-1] + h*(S0/4 *(1-alpha) + P_hum - (1-x_non_co2-x_co2)*sigma*T[-1]**4)
    if abs(T_new-T[-1]) < tolerance:
      T.append(T_new)
      break
    T.append(T_new)
  return T[-1]

  
  
if __name__=='__main__':
	co2_arr=np.arange(388*0.25, 388*4, step=int(388./4) )
	OLR_arr=[]
	for co2 in co2_arr:
		OLR_arr.append(calc_OLR(T_S=288., co2=co2))
		print("For co2=", co2," --> OLR=", OLR_arr[-1])
	T_wilson_C=calc_T(co2=co2_today, P_hum=0 )
	T_wilson_CP=calc_T(co2=co2_today, P_hum=0.034)
	print("Temperature difference for P_huum=0.034 is dT=", T_wilson_CP-T_wilson_C, " ( absolute temperatures are : ",T_wilson_CP, T_wilson_C)

#INTRO
import numpy as np
import matplotlib.pyplot as plt

import time
start_time=time.time()

import os

from EBM1D_main import area_mean
from EBM1D_main import get_D
from EBM0D import A_cald, B_cald, xi
from EBM1D_main import get_OLR
from EBM1D_main import get_S
from EBM1D_main import calc_diff
from EBM1D_main import get_albedo
from EBM1D_main import get_iceline
#from EBM1D_main import run
from EBM1D_main import get_P_hum_distribution

print("########################")
print("####  Timeseries   #####")
print("########################")

# parameters
dt=3600. *1.
nr_years=15
N=int(nr_years*3600*24*365 /dt)
C=4e7

a0=0.3
ai=0.6

#### grid
Nlats=180
lats_star=np.linspace(-90, 90, num=Nlats+1)
x_star=np.sin(lats_star*np.pi/(Nlats))
lats=(lats_star[:-1]+ lats_star[1:])/2
x=np.sin(lats*np.pi/180.)
dx=np.diff(x)
dx_star=np.diff(x_star)

T0=280-10*(3*x**2-1)
S=get_S()


#####	Equilibrium MAIN RUN FUNCTION	####
def run(T_initial, OLR_type, T_f, a0, ai, D_type, P_type, P_mean):
	Tf=T_f+273.15	# deg to C
	T=T_initial
	D=get_D(Tf,D_type, OLR_type)
	P_hum=get_P_hum_distribution(P_mean, P_type)
	T_saver=T
	for n in range(0,N):
		h = calc_diff(T, D)
		alpha=get_albedo(T, Tf, a0, ai)
		iceline=get_iceline(alpha, ai)
		dTdt=S*(1-alpha) - get_OLR(T, OLR_type)  + h + P_hum
		T=T+ dt/C * dTdt
		if (n*dt/365/24/3600) == nr_years-1:
			if area_mean(T-T_saver, lats) >0.01:
				print("After "+str(nr_years-1)+" years there's still dT>0.01 --> increase N!")
	return T, iceline
#####	DYNMAIC/TRANSIENT  MAIN RUN FUNCTION	####
def run_dyn(T_initial, OLR_type, T_f, a0, ai, D_type, P_type, P_mean):
	Tf=T_f+273.15	# deg to C
	T=T_initial
	D=get_D(Tf,D_type, OLR_type)
	P_hum=get_P_hum_distribution(P_mean, P_type)
	S=get_S()
	T_saver=T
	for n in range(0,int(float(N)/nr_years)):
		h = calc_diff(T, D)
		alpha=get_albedo(T, Tf, a0, ai)
		iceline=get_iceline(alpha, ai)
		dTdt=S*(1-alpha) - get_OLR(T, OLR_type)  + h + P_hum
		T=T+ dt/C * dTdt
	return T, iceline



# OLR
OLR_type="LM" 	# "B","LM" or "LC"

# ALBEDO/ICE
T_f=-10   ; Tf=T_f+273.15 		# -10 or -2
a0=0.3	# ice free albedo
ai=0.6	# ice albedo

# INITIAL CONDITION
T0=280-10*(3*x**2-1)  # initial Temperature

# TUNED D VALUES
D_type="D" 	# "D" or "sD"


# AHF DISTRIBUTION
P_type="G" 	# "G" or "M"

print("####  MAIN   #####")
N_time_series=150
N_cut_off=130
P_glob_today=0.034
P_glob=P_glob_today


# Uncomment if you only wnat to simulate ONE configuration
#~ OLR_type="LM"; T_f=-10; D_type="D"; P_type="G"; 
#~ P_glob=P_glob_today
#~ T_today,iceline_today= run(T0, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today)
#~ T_today,iceline_today= run(T_today, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today) # run twice to make sure that I'm in equilibrium
#~ print('initial run is finished')
#~ T_10x=T_today
#~ print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_today)+"): <T>_theta= "+'%.4f' % area_mean(T_today, lats) +", icelines:"+str(iceline_today[0])+", "+str(iceline_today[1]))

for equ_or_dyn in ["dyn"]:   #["equ", "dyn"]
	for OLR_type in ["LM"]:
		for T_f in [-10]:
			for D_type in ["D", "sD"]:
				Tf=T_f+273.15
				T_today,iceline_today= run(T0, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today)
				T_today,iceline_today= run(T_today, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today) # run twice to make sure that I'm in equilibrium
				T_10x=T_today
				folder=OLR_type+D_type+str(T_f)+P_type+"/"
				directory = os.path.dirname(folder)
				try:
					os.stat(directory)
				except:
					os.mkdir(directory)  
				print(folder, equ_or_dyn)
				t=open(folder+"T_timeseries_linear"+equ_or_dyn+".txt", 'w')
				t.write("# MODEL "+OLR_type+D_type+str(T_f)+P_type+"\n")
				for n in range(0,N_time_series):
					print("Year ", n)
					if n>N_cut_off:
						P_glob=0
					else:
						P_glob=P_glob*1.02
					
					if equ_or_dyn=="dyn":
						T_10x,iceline_10x= run_dyn(T_10x, OLR_type, T_f, a0, ai, D_type, P_type, P_glob)
					elif equ_or_dyn=="equ":
						T_10x,iceline_10x= run(T_10x, OLR_type, T_f, a0, ai, D_type, P_type, P_glob)

					print("Difference <T_10x - T_today>_theta = "+'%.6f' % area_mean(T_10x-T_today, lats)," icelines:"+str(iceline_10x[0])+" "+str(iceline_10x[1]))
					t.write(str(n)+"   "+str(P_glob)+"     "+str(area_mean(T_10x-T_today, lats))+"      "+"    "+str(iceline_10x[0])+"   "+str(iceline_10x[1]))
					t.write("\n")
				t.close()
				
print(time.time() - start_time, "  seconds needed")

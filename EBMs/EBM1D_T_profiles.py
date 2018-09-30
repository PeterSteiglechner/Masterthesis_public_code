import numpy as np
import matplotlib.pyplot as plt
import os

from EBM0D import A_cald, B_cald, xi

from EBM1D_main import area_mean
from EBM1D_main import get_D
from EBM1D_main import get_OLR
from EBM1D_main import get_S
from EBM1D_main import calc_diff
from EBM1D_main import get_albedo
from EBM1D_main import get_iceline
from EBM1D_main import run
from EBM1D_main import get_P_hum_distribution

print("########################")
print("####   T_profiles   ####")
print("########################")
#IMPLEMETNATION
C=4e7	# [J/m^2/K] heat capacity
###   Spatial Grid   ###
Nlats=180 	# Nr. of latitudes for T(phi)
lats_star=np.linspace(-90, 90, num=Nlats+1)	# latidue boundaries
x_star=np.sin(lats_star*np.pi/(Nlats)) 		
lats=(lats_star[:-1]+ lats_star[1:])/2		#latitude bands
x=np.sin(lats*np.pi/180.)
dx=np.diff(x)
dx_star=np.diff(x_star)
###  Temporal Grid  ###
dt=3600. *1. # [s]  time step for integration
nr_years=15  # [years] integrated until equilibrium
N=int(nr_years*3600*24*365 /dt)	# [s] -"-

# ALBEDO/ICE
a0=0.3	# ice free albedo
ai=0.6	# ice albedo

# INITIAL CONDITION
T0=280-20*(0.5*(3*x**2-1))  # initial Temperature


# RESULTS WITH ICE

# Choose one of the options below
#OLR_type="B"; T_f=-10; D_type="D"; P_type="G"; P_glob_today=0.034; 
#OLR_type="B"; T_f=-10; D_type="sD"; P_type="G"; P_glob_today=0.034; 
OLR_type="LM"; T_f=-10; D_type="D"; P_type="G"; P_glob_today=0.034; 
#OLR_type="LM"; T_f=-10; D_type="sD"; P_type="G"; P_glob_today=0.034; 
T_today,iceline_today= run(T0, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today)
print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_today)+"): <T>_theta= "+'%.4f' % area_mean(T_today, lats) +", icelines:"+str(iceline_today[0])+", "+str(iceline_today[1]))

for OLR_type in ["LM", "B"]:
	for T_f in [-10]:
		for D_type in ["D", "sD"]:
			P_type="G"; P_glob_today=0.034;
			print(OLR_type, T_f, D_type, P_type, P_glob_today)

			T_today,iceline_today= run(T0, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today)
			print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_today)+"): <T>_theta= "+'%.4f' % area_mean(T_today, lats) +", icelines:"+str(iceline_today[0])+", "+str(iceline_today[1]))
			P_glob_10x=0.34; 
			T_10x,iceline_10x= run(T_today, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_10x)
			print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_10x)+"): <T>_theta= "+'%.4f' % area_mean(T_10x, lats) +", icelines:"+str(iceline_10x[0])+", "+str(iceline_10x[1]))
			print("Difference <T_10x - T_today>_theta = "+'%.6f' % area_mean(T_10x-T_today, lats))
			
			# SAVE the temperature profiles
			folder=OLR_type+D_type+str(T_f)+P_type+"/"
			directory = os.path.dirname(folder)
			try:
				os.stat(directory)
			except:
				os.mkdir(directory)  
			for T in [T_10x, T_today]:
				if T[0]==T_10x[0]: 
					P=0.34
				else:
					P=0.034
				f=open(folder+OLR_type+D_type+str(T_f)+P_type+"_T_over_theta"+str(P)+".txt", 'w')
				f.write("#Model"+OLR_type+D_type+str(T_f)+P_type+", "+str(P)+'\n')
				print("###########    P_glob = "+str(P)+ "   ##############")
				# ### SAVE_DATA:
				f.write("#T_arr_"+str(T_f)+", Here: row=1 year; column= 1 lat"+'\n')
				for T_lat in T:
					f.write(str(T_lat))
					f.write("   ")
					f.write("\n")
				f.close()
	
	


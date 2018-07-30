#INTRO
import numpy as np
import matplotlib.pyplot as plt





#SOLAR
def get_S():
	# solar insolation (gridpoints x) in W/m^2
	return 1361./4*(1-0.477*0.5*(3*x**2-1))  






#MERIDIONAL TRANSPORT
def get_D(Tf, D_type, OLR_type):
	# meridional transport constant
	# will be overwritten by a more complex function later
	# For now, choose globally assigned value e.g.:
	D=global_D
	return D




def calc_diff(T, D):
	dTdx=np.diff(T)/dx
	# boundary condition: South/North Pole dTdx=0
	dTdx=np.insert(dTdx,0,0.)
	dTdx=np.append(dTdx,0.)
	h=D*np.diff(dTdx*(1-x_star**2))/dx_star
	return h


#ALBEDO
def get_albedo(T, Tf, a0, ai):
	alpha=np.array(a0+0*T)
	ind=np.where(T<Tf)
	alpha[ind]=ai
	return alpha[:]





def get_iceline(alpha, ai):
	ice_index=np.where(alpha==ai)
	icelats=lats[ice_index[0]]
	if alpha[-1]!=ai:
		iceline_N=lats[-1] # ice free NH
	else: 
		iceline_N=icelats[np.where(lats[ice_index]>0)[0][0]] # first pos icelat
	if alpha[0]!=ai:
		iceline_S=lats[0] # ice free SH
	else:
		iceline_S=icelats[np.where(lats[ice_index]<0)[0][-1]] # last neg icelat
	return [iceline_S, iceline_N] 





#OLR ###   OLR parameters   ###
##  Boltzmann  ##
##  Linear  ##
def get_OLR(T, OLR_type, co2=400):
	if OLR_type=="B":
		sigma=5.67e-8	# Boltzmann const
		tau=0.61	# tuned atm transmissivity for OLR
		return sigma*tau*T**4
	elif OLR_type=="LM":
		# Mcguffie and Henderson-Sellers (2014)
		A_M=-388.74	# W/m^2
		B_M=2.17	# W/m^2/K
		return A_M+B_M*T
	elif OLR_type=="LC":
		# Caldeira and Kasting (1992)
		from EBM0D import A_cald, B_cald, xi
		A_C=A_cald(co2)	#-324.03	# W/m^2 for 400ppm
		B_C=B_cald(co2)	#1.94	# W/m^2/K for 400ppm
		return A_C+B_C*T



#AHF
def get_P_hum_distribution(P_glob, P_type):
	if P_type=="M":
		return P_glob+0*lats
	elif P_type=="G":
		# The parameters below are from the gaussian fit to the data in Flanner (2009)
		a=0.16180316; b= 39.67767011; c=16.63696653 
		distr= a*np.exp(-(lats-b)**2/c**2)
		fit_normalized=P_glob/area_mean(distr, lats)*distr
		return fit_normalized

	
	
	
	
	
	
	
	
	
#IMPLEMETNATION
C=4e7	# [J/m^2/K] heat capacity
#####	MAIN RUN FUNCTION	####
def run(T_initial, OLR_type, T_f, a0, ai, D_type, P_type, P_mean, co2=400):
	Tf=T_f+273.15	# deg to C
	T=T_initial
	D=get_D(Tf,D_type, OLR_type)
	P_hum=get_P_hum_distribution(P_mean, P_type)
	S=get_S()
	for n in range(0,N):
		h = calc_diff(T, D)
		alpha=get_albedo(T, Tf, a0, ai)
		iceline=get_iceline(alpha, ai)
		dTdt=S*(1-alpha) - get_OLR(T, OLR_type, co2=co2)  + h + P_hum
		T=T+ dt/C * dTdt
	return T, iceline




#IMPLEMENTATION
def area_mean(variable, lats):
	weights=np.cos(lats*np.pi/180)
	return sum(weights*variable)/sum(weights)






#IMPLEMENTATION
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







#RESULTS WITHOUT ICE
T0_noIce=np.array([288 for phi in x])
global_D=0.55
T,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "G", 0.00)
T_P,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "G", 0.340)
print("Run without ice-albedo-feedback: dT_mean="+'%.5f' % (area_mean(T_P-T, lats)))

##RESULTS WITHOUT ICE
#global_D=0.55
#T,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "M", 0.00)
#T_P,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "M", 0.340)
#print("Run without ice-albedo-feedback: T_mean="+'%.5f' % (area_mean(T_P-T, lats)))

##RESULTS WITHOUT ICE
#global_D=0
#T,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "M", 0.00)
#T_P,XX= run(T0_noIce, "B", 0, 0.3, 0.3, "", "M", 0.340)
#print("Run without ice-albedo-feedback: T_mean="+'%.5f' % (area_mean(T_P-T, lats)))


#TUNING
####    Choose Model Configuration    ####
# OLR
OLR_type="B" 	# "B","LM" or "LC"

# ALBEDO/ICE
T_f=-10 		# -10 or -2
a0=0.3	# ice free albedo
ai=0.6	# ice albedo

# INITIAL CONDITION
T0=280-10*(3*x**2-1)  # initial Temperature

# TUNED D VALUES for all (OLR_type, D-type+T_f)-combinations
D_type="D" 	# "D" or "sD"
def get_D(Tf, D_type, OLR_type): # overwriting the previous simple function.
	D_B={"D-2":0.69 , "sD-2":0.665 , "D-10":0.316 , "sD-10":0.281 }
	D_LM={"D-2":0.7985, "sD-2": 0.791, "D-10":0.418 , "sD-10":0.395 }
	D_LC={"D-2":0.792, "sD-2":0.787 , "D-10":0.433 , "sD-10":0.415 }
	if OLR_type=="B":
		return D_B[D_type+str(T_f)]
	if OLR_type=="LM":
		return D_LM[D_type+str(T_f)]
	if OLR_type=="LC":
		return D_LC[D_type+str(T_f)]

# AHF DISTRIBUTION
P_type="M" 	# "G" or "M"



# RESULTS WITH ICE
OLR_type="LC"; T_f=-2; D_type="sD"; P_type="M"; P_glob_today=0.034
T_today,iceline_today= run(T0, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_today)
print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_today)+"): <T>_theta= "+'%.4f' % area_mean(T_today, lats) +", icelines:"+str(iceline_today[0])+", "+str(iceline_today[1]))
P_glob_10x=0.34
T_10x,iceline_10x= run(T_today, OLR_type, T_f, a0, ai, D_type, P_type, P_glob_10x)
print(OLR_type+D_type+str(T_f)+" ("+P_type+", P_glob="+str(P_glob_10x)+"): <T>_theta= "+'%.4f' % area_mean(T_10x, lats) +", icelines:"+str(iceline_10x[0])+", "+str(iceline_10x[1]))
print("Difference <T_10x - T_today>_theta = "+'%.6f' % area_mean(T_10x-T_today, lats))





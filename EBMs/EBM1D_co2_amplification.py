import time
start_time=time.time()
import numpy as np
import matplotlib.pyplot as plt
import os

from EBM1D_main import area_mean
from EBM1D_main import get_D
from EBM0D import A_cald, B_cald, xi
from EBM1D_main import get_OLR
from EBM1D_main import get_S
from EBM1D_main import calc_diff
from EBM1D_main import get_albedo
from EBM1D_main import get_iceline
from EBM1D_main import run
from EBM1D_main import get_P_hum_distribution

print("########################")
print("## CO2_amplification ###")
print("########################")

##### parameters######
dt=3600. *1.
nr_years=15
N=int(nr_years*3600*24*365 /dt)
C=4e7
sigma=5.67e-8
tau=0.61
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
Q=get_S()

def plotter(name, y, icelat):
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax2=ax.twinx()
	ax.plot(lats, y)
	ax.plot([icelat[1], icelat[1]], [min(y), max(y)], '--b')
	ax.plot([icelat[0], icelat[0]], [min(y), max(y)], '--b')
	plt.savefig(folder+OLR+D_type+str(T_f)+P_type+name+".png")
	return
	
#####################     MAIN      ###########################
D_type="sD" # "sD" "D"
OLR="LC" #"LM""LC""B"
T_f=-10; Tf=T_f+273.15

P_type="G" 

CO2_i=388
CO2_f= 582 #776

folder=OLR+D_type+str(T_f)+P_type+"/"
directory = os.path.dirname(folder)
try:
	os.stat(directory)
except:
	os.mkdir(directory)  

D=get_D(Tf,D_type, OLR)
print("Tf:  ", Tf,", D: ", D, ", OLR:" , OLR, ", P_type:", P_type)
#################  Spin-Up run   ##############################
print("#### Reference Run ####")
T, icelines= run(T0, OLR, T_f, a0, ai, D_type, P_type, 0 , co2=400)
print("glob_mean T", area_mean(T, lats), "iceline", icelines)
plotter("today_T_profile", T, icelines)
f=open(folder+OLR+D_type+str(T_f)+P_type+"_today_equ.txt", 'w')
f.write("T, iceline"+"\n")
f.write(str(area_mean(T, lats))+"    "+str(icelines[0])+ "   "+ str(icelines[1]))
f.close()

#################  Reference run with Pglob =0.034   ##############################
print("#### P=0.034 ####")
T_ref, icelines_ref= run(T, OLR, T_f, a0, ai, D_type, P_type, 0.034, co2=400)
print("glob_mean T", area_mean(T_ref, lats), "icelat", icelines_ref, ", DIFF: ", area_mean(T_ref-T, lats))
plotter("diff_P0034C388", T_ref-T, icelines_ref)

########## P_hum=0.34 and 3.4 and CO2=388 and 582 (all combinations  #############
print("### Check specific P-hum values in equilibrium. startT=T_P")
f=open(folder+OLR+D_type+str(T_f)+P_type+"_T_diffs_Px34CO2.txt", 'w')
f.write("#Area_mean difference of T_ref to T_Px34(P=x*0.034) using setting"+OLR+D_type+str(T_f)+P_type+'\n')
for x34 in [0.034, 0.34]:
	for co2 in [CO2_i, CO2_f]:
		print("###########    P_glob = "+str(x34)+ ", CO2="+str(co2)+"   ##############")
		T_x34=T_ref 
		T_x34, icelines_x34 = run(T_x34, OLR, T_f, a0, ai, D_type, P_type, x34 , co2=co2)
		f.write("P_hum="+str(x34)+", CO2="+str(co2)+"\n")
		f.write(str(x34)+"   "+str(area_mean(T_x34-T_ref, lats))+"    "+str(icelines_ref[0])+","+str(icelines_ref[1])+"    "+str(icelines_x34[0])+","+str(icelines_x34[1])+'\n')
		print("P", x34, "C", co2, ", DIFF: ", area_mean(T_x34-T_ref, lats), " icelat: ", icelines_x34)
		plotter(str(co2)+str(x34), T_x34, icelines_x34)
		plotter("diff"+str(co2)+str(x34), T_x34-T_ref,icelines_x34)
f.close()

print(time.time() - start_time, "  seconds needed")


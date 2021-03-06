import time
start_time=time.time()
import numpy as np
import matplotlib.pyplot as plt
import os

import EBM1Dmodel as model
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
print("## Table_tenfold_Phum ##")
print("########################")

D_types=["D", "sD"] # "sD
OLRs=["LM", "LC", "B"] #"LM""LC"
T_fs=[-2, -10]
P_types=["M", "G"] 

# parameters
model.dt=1./4*3600. *1.
model.nr_years=15
model.N=int(model.nr_years*3600*24*365 /model.dt)
model.C=4e7

a0=0.3
ai=0.6

#### grid
model.Nlats=360
model.lats_star=np.linspace(-90, 90, num=model.Nlats+1)
model.x_star=np.sin(model.lats_star*np.pi/180)
model.lats=(model.lats_star[:-1]+ model.lats_star[1:])/2
model.x=np.sin(model.lats*np.pi/180.)
model.dx=np.diff(model.x)
model.dx_star=np.diff(model.x_star)

T0=280-10*(3*model.x**2-1)
Q=get_S()  #solar insolation (lats) in W/m^2


for D_type in D_types:
	for OLR in OLRs:
		for T_f in T_fs:
			for P_type in P_types:
				folder=OLR+D_type+str(T_f)+P_type+"_doubleRes/"
				directory = os.path.dirname(folder)
				try:
					os.stat(directory)
				except:
					os.mkdir(directory)  
				print(folder)

				#####################     MAIN      ###########################
				Tf=T_f+273.15
				D=get_D(Tf,D_type, OLR)
				print("Tf:  ", Tf,", D: ", D, ", OLR:" , OLR, ", P_type:", P_type)
				#################  Reference run with Pglob =0.034   ##############################
				print("#### P=00.034 ####")
				T_P, icelines_P = run(T0, OLR, T_f, a0, ai, D_type, P_type, 0.034 , co2=400)
				print("glob_mean T", area_mean(T_P, model.lats), "icelines", icelines_P)
				### SAVE_DATA:
				f=open(folder+OLR+D_type+str(T_f)+P_type+"0.034.txt", 'w')
				f.write("T, iceline"+"\n")
				f.write(str(area_mean(T_P, model.lats))+"    "+str(icelines_P[0])+ "   "+ str(icelines_P[1]))
				f.write("\n")
				f.close()

				############### Perturbed RUN   ##########################
				########## P_hum=0.34 and 3.4 and integrate to equilibrium  #############
				print("### Check specific P-hum values in equilibrium. startT=T_P")
				for x34 in [0.34]:
					f=open(folder+OLR+D_type+str(T_f)+P_type+"_T_diffs_P"+str(x34)+".txt", 'w')
					f.write("#Area_mean difference of T_P to T_Px34(P=x*0.034) using setting"+OLR+D_type+str(T_f)+P_type+'\n')
					print("###########    P_glob = "+str(x34)+ "   ##############")
					T_x34, icelines_x34 = run(T_P, OLR, T_f, a0, ai, D_type, P_type, x34 , co2=400)
					f.write(str(x34)+"   "+str(area_mean(T_x34-T_P, model.lats))+",     "+str(icelines_P[0])+","+str(icelines_P[1])+",   "+str(icelines_x34[0])+","+str(icelines_x34[1])+'\n')
				f.close()


print("###################################################")
print("## Transfer to format, that can be used in latex ##")
print("###################################################")
import StringIO

growth_rate=1.02
icelat_S_i={}
icelat_S_f={}
icelat_N_i={}
icelat_N_f={}
dT={}
row=0

for D_type in D_types:
    for OLR in OLRs:
        for T_f in T_fs:
            for P_type in P_types:
                folder=OLR+D_type+str(T_f)+P_type+"_doubleRes/"
                filename=folder+OLR+D_type+str(T_f)+P_type+"_T_diffs_P0.34.txt"
                s = open(filename).read().replace(',','   ')
                data=np.loadtxt(StringIO.StringIO(s))
                name=OLR+D_type+str(T_f)+P_type
                print(name)
                icelat_S_i[name]=data[2]
                icelat_S_f[name]=data[4]
                icelat_N_i[name]=data[3]
                icelat_N_f[name]=data[5]
                dT[name]=data[1]
tot={}
tot["icelat_S_i"]=icelat_S_i
tot["icelat_S_f"]=icelat_S_f
tot["icelat_N_i"]=icelat_N_i
tot["icelat_N_f"]=icelat_N_f
tot["dT"]=dT
def print_setting(tot,OLR ,D_type, T_f,P_type, verbose=True):
    name=OLR+D_type+str(T_f)+P_type
    string=r"$\mathbf{"+'%.3f' % (tot["dT"][name])+"\, \degC}$  &  $["+ str(tot["icelat_S_i"][name])+" , "+str(tot["icelat_N_i"][name])+r"]$ \arrow $["+ str(tot["icelat_S_f"][name])+" , "+str(tot["icelat_N_f"][name])+"]$"
    if verbose:
        print("Model "+name+": Changes due to P_hum=0.34")
        print("dT: ", tot["dT"][name])
        print("icelat:   ["+ str(tot["icelat_S_i"][name])+" , "+str(tot["icelat_N_i"][name])+"] \\arrow ["+ str(tot["icelat_S_f"][name])+" , "+str(tot["icelat_N_f"][name])+"]")
        print("")
        print(string)
    return string
#print_setting(tot, "B", "sD", -2, "M", verbose=False)                
def print_full_table(tot,OLR ,D_type, T_f,P_type, verbose=True):
    if P_type=="G":
        _P="G\\textsubscript{F}"
    else:
        _P=P_type
    pre_string="("+OLR+D_type+"$_{"+str(T_f)+"}$("+_P+"))&"
    #(LMsD$_{-10}$(G\textsubscript{F})&
    string=pre_string+print_setting(tot,OLR ,D_type, T_f,P_type, verbose=False)+"\\"+"\\"
    return string
#print_full_table(tot, "B", "sD", -2, "M", verbose=False)

full_table=open("EBM1D_table_tenfold_doubleRes.txt", 'w')
full_table.write("\\textbf{Model Configuration} & $\\mathbf\{<\\Delta T>_\\theta}$ & \\textbf{initial} \\arrow \\textbf{final iceline [South, North]} \\\\ \\hline "+"\n")
for OLR in OLRs:
    for D_type in D_types:
        for T_f in T_fs:
            for P_type in P_types:
                full_table.write(print_full_table(tot, OLR, D_type, T_f,P_type, verbose=False))
                full_table.write("\n")
full_table.close()

print(time.time() - start_time, "  seconds needed")


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ncep_url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/"
ncep_Ts = xr.open_dataset( ncep_url + "surface_gauss/skt.sfc.mon.1981-2010.ltm.nc", decode_times=False)
ncep_ulwrf = xr.open_dataset( ncep_url + "other_gauss/ulwrf.ntat.mon.1981-2010.ltm.nc", decode_times=False)
olr_lonAvg=ncep_ulwrf.mean(dim='lon')
Ts_lonAvg=ncep_Ts.mean(dim='lon')

fs=27
font={'family' : 'serif', 'size' : fs}		#define font for legend
plt.rc('font', **font)	
plt.rc('text', usetex=True)

import sys
sys.path.append("../../../EBM_equilibrium/")
from EBM0D_Wilson_and_GeaBanacloche2012 import *
#sys.path.append("../../../Woche9/")
#import wilson_class
#a=wilson_class.wilson()
#T_arr=np.linspace(260,305, num=25)
#OLR=a.calc_OLR(T_S=T_arr, co2=400)
#T_mit=a.calc_T(co2=400, P_hum=0.034)
#T_ohne=a.calc_T(co2=400, P_hum=0)
#print("Wilson delta T=", T_mit-T_ohne)
T_mit=calc_T(co2=400, P_hum=0.034)
T_ohne=calc_T(co2=400, P_hum=0)
print("Wilson delta T=", T_mit-T_ohne)
T_arr=np.linspace(210,310, num=40)
OLR=[]
for i in T_arr:
	OLR.append(calc_OLR(T_S=i, co2=400, co2_today=400))
print(OLR)	

fig=plt.figure(figsize=(16,9))
ax=fig.add_subplot(111)
ax.scatter(Ts_lonAvg.skt+273.15, olr_lonAvg.ulwrf, alpha=0.4,label=r'Monthly mean data for'+'\n'+'OLR and $T_S$ 1981-2010')
T=np.linspace(210, 310, num=2)
#ax.plot(T, 210-273.15*2 +2.*T, '-g', lw=3, label=r'(Rose, 2017)')
ax.plot(T, -388.74+2.17*T, '-r', lw=3, label=r'(Mcguffie and Hendersson-Sellers, 2014)')
#import sys
#sys.path.append("../../../Woche9/")
#import caldeira_class as cald
#a=cald.caldeira()
#A_cald=a.A_cald(400)
#B_cald=a.B_cald(400)
sys.path.append("../../../EBM_equilibrium/")
from EBM0D import *
A_cald= A_cald(400)
B_cald= B_cald(400)

ax.plot(T, A_cald+B_cald*T, '-k', lw=3, label=r'(Caldeira and Kasting, 1992)'+'\n'+' for $CO_2=400\, ppm$')
m=0.5 #m = atmospheric transmission factor
Sellers_OLR=5.67e-8 * T**4 *(1-m*np.tanh(19*T**6*1e-16))
ax.plot(T, Sellers_OLR, '-', color='m', lw=3, label=r'(Sellers, 1969)'+'atm. transm. $m=$'+str(m))
umrechnungs_faktor=4184./(24*3600*30)*10000 # von kcal/cm^2 to W/m^2
a=14; B=0.14; a_1=3; B_1=0.1 # In kcal/cm^2 bzw. kcal/cm^2/degC
a=umrechnungs_faktor*(a - 273.15*B)
a_1=umrechnungs_faktor*(a_1-273.15*B_1)
B=umrechnungs_faktor*B
B_1=umrechnungs_faktor*B_1
n=0.5 #n=cloud_fraction
ax.plot(T, a+B*T-n*(a_1+B_1*T), '-c', lw=3, label=r'(Budyko, 1969), cloud fraction $n=$'+str(n))
ax.plot(T_arr, OLR, '--', color='orange', lw=5, label=r'(Wilson and Gea-Banacloche, 2012) '+'\n'+'non-linear model, $CO_2=400\, ppm$')
ax.plot(T_arr, 0.61*5.67e-8*T_arr**4, '--', color='green', lw=5, label=r'Boltzmann Radiation '+r'$\tau = 0.61$')

ax.legend(fontsize=21)
ax.set_ylabel(r"$OLR \ [W/m^2]$")
ax.set_xlabel(r"$T_S \ [K]$")
ax.grid()
ax.set_xlim(210, 305)
ax.set_ylim(100, 300)
#ax.annotate(r"More (badly fitted) $OLR(T_S)$-data for $T_S<260\, K$", 
#            xy=(260,170), xytext=(275, 170), arrowprops={'facecolor':'blue', 'alpha':0.4})
plt.tight_layout()
plt.savefig("analyse_linear_olrs.png")
plt.show()




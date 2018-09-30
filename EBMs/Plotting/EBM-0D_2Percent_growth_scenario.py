#9.2.18 Peter Steiglechner
#SAME PROGRAMM AS RADIATIVE_BALANCE_BERG but new


import numpy as np
from matplotlib import pyplot as plt
fs=22						#fontsize for labels, legend, ...
font={'family' : 'serif', 'size' : fs}		#define font for legend
plt.rc('font', **font)	

sig=5.67e-8 #W/m**2 /K**4
alpha=0.3
S0=1361. #W/m**2
tau=0.61

T_S_no_atm=(S0/4 *(1-alpha) / (sig))**(1./4)


print("Temperature without atmosphere is: ",  '%.3f' % T_S_no_atm, " bzw. ",  '%.3f' % (T_S_no_atm-273.15))

T_S_w_atm=((S0/4 *(1.-alpha))/ (sig* tau))**(1./4) 
A=-388.74
B=2.17
T_S_w_atm_linear=(S0/4 *(1-alpha)-A) / B

glob_mean_T_preind=np.copy(T_S_w_atm)	#for further use

#INTRODUCE P_hum:
P_hum_now=0.036

# RAD Difference: 
#print("radiative Difference: ", '%.3f' % P_hum, "W/m^2")
def T_AHF(P_hum):
	return ((S0/4 *(1.-alpha)+P_hum)/ (sig* tau))**(1./4)
def T_AHF_lin(P_hum):
	return ((S0/4 *(1.-alpha)+P_hum-A)/ B)

#T_S_AHF=T_AHF(P_hum)

#############################################
# Here, i plot foor different growth scenarios (2% or 2.3% Power increase per yr) the resulting Temp growth due to AHF. 
# P_hum at present needs to be inserted as total power consumption [W] that dissipates as heat divided by 4pi r^2
#This assumes that P_hum comes from entirely external energy sources (fossil fuels)
def plot_time_evolution(P_hum_now, growth_rate, ax_T, ax_P, line, text):
	N=200
	l=np.empty([N, 4])			
	l[:,0]=np.linspace(2014,2014+N, num=N)		#years
	P_hum=np.copy(P_hum_now)
	for i in range (0,N):
		P_hum=P_hum*(1+growth_rate/100)	 #increase of 2.9% per year
		l[i,2]=np.copy(P_hum)
		l[i,1]=T_AHF(P_hum)
		l[i,3]=T_AHF_lin(P_hum)

	a, =ax_T.plot(l[:,0], l[:,1]-glob_mean_T_preind, label="$\Delta T_S$ for "+str(growth_rate)+"%-growth per year", lw=5, color="blue", ls=line)
	b, =ax_P.plot(l[:,0], l[:,2], label=r"$P_{\rm hum}$ for "+str(growth_rate)+"%-growth per year", lw=2, color="red", ls=line)
	#c, =ax_T.plot(l[:,0], l[:,3], label="Temperature (linear) for "+str(growth_rate)+"%/a -- "+text, lw=1, color="red", ls=line)

	ax_P.set_ylim(bottom=0)

	ax_T.set_xlabel("years", fontsize=fs)
	ax_T.tick_params(axis='y', colors='blue')
	ax_P.tick_params(axis='y', colors='red')
	ax_T.set_ylabel(r"$ \overline{\Delta T_S} \,  [{\rm K}]$", fontsize=fs, color="blue")
	ax_P.set_ylabel(r"$P_{\rm hum} \, [{\rm W/m^2}]$", fontsize=fs, color="red")
	return a, b
	

text='fossil_fuels'
fig= plt.figure(figsize=(16,9))
ax1=fig.add_subplot(111)
ax2=ax1.twinx()
a_29, b_29= plot_time_evolution(P_hum_now, 2.3, ax1, ax2, '--', text)
a_20, b_20= plot_time_evolution(P_hum_now, 2.0, ax1, ax2, '-', text)

data1=np.loadtxt("world_only_Energy_use.txt")
years1=data1[0]
energy_use=data1[1][:]
data2=np.loadtxt("world_only_Pop.txt")
years2=data2[0][:-2]
pop=data2[1][:-2]
if len(years2)!=len(years1):
    print("not the same length")
umrechnung=11.63 #kWh/a/Kopf   pro kg oel
tot_E=umrechnung*energy_use*pop*1e3
P_hum_past=tot_E*1./(365*24)/(4*np.pi*6371000**2)
ax2.plot(years1, P_hum_past, 'or', ms=10, label=r"historic $P_{\rm hum}$-data"+" from The World Bank \n (`Energy Use per Capita' and `total Population')")
import scipy.optimize as opt
def f(y, a, b):
    return (1.00+a)**(y-y[0]) *b 
popt, pcov=opt.curve_fit(f, years1, P_hum_past, bounds=([0.00, 0], [1, np.inf]))

left, bottom, width, height = [0.2, 0.25, 0.3, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
#ax3.set_yticks(ax2.get_yticks(), minor=True)
ax3.set_yticks([0.02, 0.03])#, minor=True)
ax3.yaxis.tick_right()
ax3.tick_params(axis='y', colors='red')
#ax3.set_xticks([1971, 1985, 2000, 2014])
#ax2.grid(b=True, which='minor')
ax3.minorticks_on()
ax3.grid(b=True, which='major', alpha=0.7)
ax3.grid(b=True, which='minor', alpha=0.3)
cross,= ax3.plot(years1, P_hum_past, 'o', ms=10, color='red', label=r"$P_{\rm hum}$:"+" historic data from The World Bank \n (`Energy Use per Capita' and `total Population')")
fit, = ax3.plot(years1, f(years1, *popt), '-g', lw=2,  label=r'exp. Fit to historic $P_{\rm hum}$ with growth rate '+'%.2f' % (popt[0]*100)+"% per year")
ax3.set_xticks([1970+k*10 for k in range(0,5)])




plots=[a_20,b_20, a_29, b_29, cross, fit]

ax1.legend(plots, [i.get_label() for i in plots], loc="best", fontsize=18)
ax1.set_title(r"EBM-0D: $T_S$ increase for exponential $P_{\rm hum}$-growth"+"\n")
ax1.set_yticks(np.arange(0, 1.2, 0.2))
ax1.set_xticks(np.arange(2015, 2015+200+1, step=40))
ax1.set_xlim(1965, 2015+200+1)
ax1.grid(b=True, axis='both')

plt.savefig("T_evo_AHF_exp_growth_"+text+".pdf")
plt.show()







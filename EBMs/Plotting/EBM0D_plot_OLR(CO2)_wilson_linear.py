'''
Peter Steiglechner 9.2.
THIS PROGRAM RESEMBLES WEEK9 /home/peter/PIK/Woche9/CO2_impact_comparison_Wilson-Caldeira.py
  
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

fs=22
font={'family' : 'serif', 'size' : fs}		#define font for legend
plt.rc('font', **font)	
#plt.rc('text', usetex=True)

sys.path.append("../../../EBM_equilibrium/")
from EBM0D_Wilson_and_GeaBanacloche2012 import *

sys.path.append("../../../EBM_equilibrium/")
from EBM0D import *
A_cald_tod= A_cald(400)
B_cald_tod= B_cald(400)
T_S=288

OLR_today_wils=calc_OLR(co2=400, T_S=T_S)
INC=1361./4*(1-0.3)
OLR_today_cald=(A_cald_tod+B_cald_tod*T_S)
print(INC, OLR_today_cald, OLR_today_wils)

co2_arr=np.arange(400*0.25, 400*4, step=int(400./4) )
OLR_cald=(A_cald(co2_arr)+B_cald(co2_arr)*T_S)
OLR_wilson=[]
for c in co2_arr:
	OLR_wilson.append(calc_OLR(co2=c, T_S=288))

doubling_cald=(A_cald(2*400)+B_cald(2*400)*T_S) - (A_cald(400)+B_cald(400)*T_S)
doubling_wilson=calc_OLR(co2=2*400, T_S=288)-calc_OLR(co2=400, T_S=288)

string=r'\underline{$\Delta OLR$ for CO$_2$-doubling ($400$ to $800\, ppm$): }'+' \n '
string+=r"\begin{itemize} \item Linear:  "+ '%.4f' % doubling_cald+r'$\, W/m^2$'#+' \n '
string+=r"\item Wilson:  "+ '%.4f' % doubling_wilson+r'$\, W/m^2$'
string+=r"\item IPCC,2013 (fixed SST): $-3.7 \pm 0.8 \, W/m^2$"
string+="\end{itemize}"
print(string)

OLR_today=5.67e-8 *0.61*288**4
print("OLR", OLR_today)
OLR_ipcc=OLR_today-3.7/np.log(2)*np.log(co2_arr/400)
OLR_ipcc_upper=OLR_today-(3.7+0.8)/np.log(2)*np.log(co2_arr/400)
OLR_ipcc_lower=OLR_today-(3.7-0.8)/np.log(2)*np.log(co2_arr/400)

fig=plt.figure(figsize=(16,9))
ax=fig.add_subplot(111)
ax.plot(co2_arr, [INC for _ in co2_arr], '--', label=r'$INC=S_0/4\cdot (1-\alpha_P)$', color='green', lw=2)
ax.plot(co2_arr, [5.67e-8 *0.61*288**4 for _ in co2_arr],'--', color='orange',lw=2,
label=r'$\tau \cdot \sigma T^4 = 0.61\cdot \sigma \cdot (288\, {\rm K})^4$')

ax.plot(co2_arr, OLR_wilson, label='(Wilson and Gea-Banacloche, 2012)', color='darkblue', lw=3)
ax.plot(co2_arr, OLR_cald, label='(Caldeira and Kasting, 1992)', color='red', lw=3)


ax.plot(co2_arr, OLR_ipcc, lw=4, label=r'(IPCC, 2013) eff. rad. forcing: '+'\n'+r'$3.7\pm 0.8 \, {\rm W/m^2}$'+' per CO$_2$-doubling', color='black')
ax.fill_between(co2_arr, OLR_ipcc_upper, OLR_ipcc_lower, alpha=0.3, color='grey' )
#Table 9.5 Page 834 / 818 in IPCC AR 5

ax.set_xlabel(r"CO$_2$-concentration in ${\rm ppm}$")
ax.set_ylabel(r"OLR in ${\rm W/m^2}$")
ax.legend()
ax.grid()
#ax.annotate(string,  
#	xy=(800,239), 
#	fontsize=20)
'''
ax.annotate(r'\textbf{WILSON:} Via $CO_2$-absorption sepctrum: '+r' $x_{CO_2}(CO_2\_conc)$'+' \n '+r'$(1-(x_{non-CO_2} + x_{CO_2})) \cdot \sigma T^4$'+' \n' +r'tuning parameter e.g. '+'$x_{non-CO_2}, \ T_{Tropopause}=217K$ ', 
xy=(co2_arr[2],OLR_wilson[2]), 
xytext=(co2_arr[2]+20,OLR_wilson[2]+4), 
arrowprops=dict(facecolor='darkblue', shrink=0.02), 
fontsize=20)
ax.annotate(r'\textbf{LINEAR:} $A(ln(CO_2))+B(ln(CO_2))\cdot T$ ', 
xy=(co2_arr[2],OLR_cald[2]), 
xytext=(co2_arr[2]-80,OLR_cald[2]-6), 
arrowprops=dict(facecolor='red', shrink=0.02), 
fontsize=20)
'''
fig.tight_layout()
plt.savefig("T_CO2_comparison_linear-wilson.png")
plt.show()


#~ string2="INC="+ '%.4f' % (1361./4*(1-0.3))+'\n'+ "OLR_caldeira="+\
#~ '%.4f' % OLR_today_cald + '\n'+"OLR_wilson="+'%.4f' % OLR_today_wils

#~ print(string2)



#################################################################
##############    Compute Temperature and compare   #############
#################################################################
#T_wilson=a.calc_T(co2=450)
#print(T_wilson)
#co2_arr=np.arange(380*0.5, 380*3, step=int(380./2) )
co2_arr=np.arange(400*0.25, 400*4, step=int(400./4) )
co2_arr_saved=np.arange(400*0.25, 400*4, step=int(400./4) )
print(np.shape(co2_arr))
#~ T_wilson=a.calc_T_arr(co2_arr)
#~ f=open("saved_T_Wilson_CO2_impact_comparison.txt", "w")
#~ for i in T_wilson:
	#~ f.write(str(i)+'\n')
#~ f.close()
print(co2_arr, co2_arr_saved, (co2_arr==co2_arr_saved).all)
if (co2_arr==co2_arr_saved).all:
	T_wilson=np.loadtxt("saved_T_Wilson_CO2_impact_comparison.txt")
else: 
	T_wilson=a.calc_T_arr(co2_arr)
print(T_wilson)
T_linear=(INC - b.A_cald(co2_arr))/b.B_cald(co2_arr)

T_S_400=288.
co2_400=400.
T_ipcc=T_S_400+3.2/np.log(2)*np.log(co2_arr/co2_400)
T_ipcc_upper=T_S_400+(3.2+1.3)/np.log(2)*np.log(co2_arr/co2_400)
T_ipcc_lower=T_S_400+(3.2-1.3)/np.log(2)*np.log(co2_arr/co2_400)

fig=plt.figure(figsize=(10,9))
ax=fig.add_subplot(111)
#ax.plot(co2_arr, T_linear, label=r'\textbf{LINEAR:} $A(ln(CO_2))+B(ln(CO_2))\cdot T$', color='red', lw=4)
ax.plot(co2_arr, T_linear, label=r'(Caldeira and Kasting, 1992)', color='red', lw=4)
#ax.plot(co2_arr, T_wilson, label=r'\textbf{WILSON:} $(1-(x_{non-CO_2} + x_{CO_2})) \cdot \sigma T^4$, '+'\n'+r'via $CO_2$-spectrum integration', color='blue', lw=4)
ax.plot(co2_arr, T_wilson, label=r'(Wilson and Gea-Banacloche, 2012)', color='blue', lw=4)
#ax.plot(co2_arr, T_ipcc, lw=4, label=r'\textbf{IPCC:} $T_S=$'+'%.1f' % T_S_400 + r" for $CO_2=$"+'%.1f' % co2_400+'\n'+
#	r'equ. cl. sens. $ECS=3.2\pm 1.3 \, K$ per $CO_2$-doubling'+' \n ' +
#	 r'$T=T_{today} + ECS \cdot ln(co_2/co_{2,today})$' , color='black')
ax.plot(co2_arr, T_ipcc, lw=4, label=r'\textbf{IPCC:} equilib. clim. sens. '+'\n'+r'$ECS=3.2\pm 1.3 \, K$'+'\n'+'per $CO_2$-doubling', color='black')

ax.fill_between(co2_arr, T_ipcc_lower, T_ipcc_upper, alpha=0.3, color='grey' )
ax.legend()
ax.set_xlabel(r"$CO_2$ concentration $[ppm]$")
ax.set_ylabel(r"Temperature $[K]$")
ax.set_title(r'$CO_2$-dependency of OLR-parametrizations')
ax.grid()
#~ ax.annotate(r'\textbf{IPCC:} equilibrium climate sensitivity $ECS=3.2\pm 1.3 \, K$ per doubling of $CO_2$'+' \n ' + r'$T=T_{today} + ECS \cdot ln(co_2/co_{2,today})$' ,
		#~ xy=(co2_arr[1],T_ipcc_upper[1]), 
		#~ xytext=(co2_arr[1]+100,T_ipcc_upper[1]-1), 
		#~ arrowprops=dict(facecolor='black', shrink=0.02), 
		#~ fontsize=18)
#~ ax.annotate(r'\textbf{WILSON:} $(1-(x_{non-CO_2} + x_{CO_2})) \cdot \sigma T^4$'+' \n' +r'via $CO_2$-spectrum integration', 
		#~ xy=(co2_arr[5],T_wilson[5]), 
		#~ xytext=(co2_arr[5]+200,T_wilson[5]-1.7), 
		#~ arrowprops=dict(facecolor='blue', shrink=0.02), 
		#~ fontsize=18)
#~ ax.annotate(r'\textbf{LINEAR:} $A(ln(CO_2))+B(ln(CO_2))\cdot T$', 
		#~ xy=(co2_arr[1],T_linear[1]), 
		#~ xytext=(co2_arr[1]-70,T_linear[1]+3.0), 
		#~ arrowprops=dict(facecolor='red', shrink=0.02), 
		#~ fontsize=18)
plt.tight_layout()
plt.savefig("CO2_comparison_linear-wilson.png")
plt.show()


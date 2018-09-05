import numpy as np
data_co2=np.loadtxt("co2_past.dat")
data_co2_fut=np.loadtxt("SSP5baseline_CO2.txt")
years_past=data_co2[:,0]
years_ssp=data_co2_fut[:,0]
co2_past=data_co2[:,1]
co2_fut=data_co2_fut[:,1]

growth_rate_per_year=(co2_fut[1:]-co2_fut[:-1])/np.diff(years_ssp)

patched_co2=[co2_past[-1]]
patched_years=[years_past[-1]]
for i in range(0,len(years_ssp[1:-1])):
    for m in range(int(years_ssp[1+i])+1, int(years_ssp[1+i+1])+1):
        last_co2=patched_co2[-1]
        patched_co2.append(last_co2+growth_rate_per_year[1+i])
        patched_years.append(m+0.5) # historic data is defined at year=YYYY.5
f=open("co2_ssp5.dat", 'w')
for n,i in enumerate(years_past):
    f.write('        '+str(years_past[n])+"        "+str(co2_past[n])+"\n")
for i in range(1,len(full_years)):
    f.write('        '+str(full_years[i])+"        "+str(full_co2[i])+"\n")
f.close()

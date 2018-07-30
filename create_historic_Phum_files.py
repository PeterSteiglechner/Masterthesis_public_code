# Peter Steiglechner, 28.2.18
# This program proudces a netcdf file or multiple txt-files including the P_hum-Matrix in each year from 1800 up to 2010
# The data is taken from the WOrld Bank (see below)
# ACTION: Adjust the path for the Energy data from the WorldBank (line 36) 
#      and (in order to get the climber grid) the path to a random snapshot of the atmosphere of climber (line 65)


from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import os

# ## Data from the World Bank
# 
# http://databank.worldbank.org/data/reports.aspx?source=2&series=EG.USE.PCAP.KG.OE&country=EAS#
# 
# The data for Energy per Capita (kg Oil/a) and population for several country aggregates is saved in an excel file/txtfile:
#   E_consumption_population_world_final_filledEUWithExtrapolatedEpC.xlsx
#   E_consumption_area_fromExcel.txt
# Via a conversion factor (kg Oil in kWh) this is converted to a total energy used in each world region

regions=['World','Sub-Saharan Africa', 'South Asia', 'North America', 'Middle East & North Africa',
         'Latin America & Caribbean', 'Europe & Central Asia','East Asia & Pacific without Australia', 'Australia']
r=['w', 'ssa', 'sa', 'na', 'mena', 'lac', 'eca', 'eap-aus', 'aus']

def get_data_of_year(array, y, years):
    h=np.empty([len(r)])
    ind=(np.where(years==y)[0][0])
    for num,i in enumerate(r):
        h[num]=array[i][ind]
    return h

def get_Energy():
    data_E=np.loadtxt(E_consumption_area_fromExcel.txt")
    years=data_E[0,:]
    EpC={}
    Pop={}
    for i in range(0,len(r)):
        Pop[r[i]]= data_E[1+2*i]
        EpC[r[i]] = data_E[2+2*i,:] #in kg oil/capita/yr
    E_tot={} 
    umrechnungs_faktor=11.63 #kWh pro kg oil #SIEHE GOOGLE
    for i in r:
        E_tot[i]=Pop[i]*EpC[i]*umrechnungs_faktor
        if i=='eap-aus':
            E_tot['eap-aus']=umrechnungs_faktor*(Pop['eap-aus']*EpC['eap-aus'] - Pop['aus']*EpC['aus'])
    perc={}
    for num, i in enumerate(r[1:]):
        perc[i]=E_tot[i]/E_tot['w']
    tot_E_sum=np.empty([len(years)])
    for num,y in enumerate(years):
        Ey=get_data_of_year(E_tot, y, years)
        tot_E_sum[num]= np.sum(Ey[1:])
    ##### Because the area_mean(PMap) is a bit smaller than E_tot_glob/Area_glob due to the 2% we forgot in the energy data
    ##### Therefore we simply divide by 0.98 and get to sligthly higher values in all regions but in total we get E_tot_glob/Area_glob
    renormalization_factor= (tot_E_sum/E_tot['w'])
    return years, perc, E_tot, tot_E_sum, renormalization_factor

# ## Data from a snapshot of a climber 3 a run in order to get land fraction
def get_climber_grid():
    data_grid=xr.open_dataset("snapshots_potsdam2.055623.01.01.dta.nc", decode_times=False)
    return data_grid

def area_lon(x):
    '''
    x should have dimensions: lat,lon and ones where the area is calculated or just the land fraction
    gives back the SIZE OF THE AREA, where x has some land fraction. 
    (for x=landfract this gives the total land fraction area =0.29)
    '''
    dlon=(x.xt_i[1]-x.xt_i[0])*np.pi/180
    dlat=(x.yt_j[1]-x.yt_j[0])*np.pi/180
    ind=np.where(x!=0) #gives the indices (max = lats*lon = 24*16) where x is not 0
    lats=ind[0]
    lons=ind[1]
    area=0
    for i in range(0,len(lats)):
        lat_ind=lats[i]
        lon_ind=lons[i]
        area+=x[lat_ind,lon_ind] * np.cos(np.pi/180*x.yt_j[lat_ind])*dlon*dlat ## fraction * cos(theta) dtheta dphi
    return area
#area_lon(data_grid.frlnd[0]*get_lat_lon_agg('na'))

def get_lat_lon_agg(aggregate, climber_grid):
    '''
    #Asign each world region to a square (and some additional) cells of CLIMBE
    #Note: this is quite random. Some visualization of that is in check_ncfile one folder up
    #Note: I have removed Russia, North Canada and Australia
    '''
    array=np.zeros_like(climber_grid)
    if aggregate=='sa': # 'South Asia':
        lat_low=13; lat_high=16
        lon_low=3; lon_high=3
        #additional=[[16,2]]
        additional=[]
    if aggregate=='lac':#'Latin America and Caribbean':
        lat_low=4; lat_high=15
        lon_low=10; lon_high=14
        additional=[]
    if aggregate=='ssa':#'Sub-Saharan Africa':
        lat_low=7; lat_high=14
        lon_low=0; lon_high=1
        additional=[[13,15], [14,15], [12, 15]]
    if aggregate=='mena':#'Middle East and North Africa
        lat_low=15; lat_high=16
        lon_low=0; lon_high=2
        additional=[[16,15],[15,15], [14,2]]	
    if aggregate=='eca':#'Europe and Central Asia
        lat_low=17; lat_high=20
        lon_low=0; lon_high=7
        additional=[[17,15],[18,15], [19,15], [20,15]]    
        ########################################### CUT OFF RUSSIA
        lon_high=2
    if aggregate=='eap-aus':#'East Asia Pacific
        lat_low=6; lat_high=17
        lon_low=4; lon_high=6
        additional=[[17,3]]  
        ######################## CUT OFF AUSTRALIA
        lat_low=10
    if aggregate=='aus':#Australia #### NOTE: REMOVE ENERGY OF AUSTRALIA FROM EAP...
        lat_low=6; lat_high=9; lon_low=4; lon_high=7
        additional=[]
    if aggregate=='na':#'NorthAMerica
        lat_low=16; lat_high=21
        lon_low=10; lon_high=12
        additional=[]    
        ########################################### CUT OFF North CANADA
        lat_high=18
    for lat in range(lat_low, lat_high+1):
        for lon in range(lon_low, lon_high+1):
            array[lat,lon]=1
    for p in additional:
        array[p[0], p[1]]=1
    return array

def get_P_map_year(year, years, data_grid, E_tot, perc, renormalization_factor):
    '''
    ## Function to create the Map of P. 
    #This takes the E_tot of the world, the percentages of each world region of the E_tot, the land (-glacier) fraction in the grid
    #For one year this calculates for each world region the totE produced there, divides this by the size of the region and 
    #divides by the forgotten islands factor (=0.979 which makes all P's a bit bigger in order to match to the global value) 
    #and hence results in a Map of P_hum
    '''
    year_ind=np.where(years==year)[0][0]
    tot=E_tot['w'][year_ind] *1000/(365*24)   #to get from kWh to W/m^2 
    P_map=np.zeros_like(data_grid.frlnd[0])
    for num, reg in enumerate(r[1:]):
        percent=perc[reg][year_ind]
        region_array=get_lat_lon_agg(reg, data_grid.frlnd[0])
        current_area=region_array*(data_grid.frlnd[0]-data_grid.frlnd[0]*data_grid.frglc[0])
        area_size=area_lon(current_area)*6371000**2
        P_map+=percent*tot*np.array(region_array) *1./ area_size *1./renormalization_factor[year_ind]  
    return P_map
    


### CREATE FULL ARRAY over the relevant time period
##### Before 1971: Assume 2 % growth until 1971 in the same distribution as in 1971
if __name__ == '__main__':
    years, perc, E_tot, tot_E_sum, renormalization_factor= get_Energy()
    
    data_grid=get_climber_grid()
    
    start_year=1810
    end_year=2010
    fullyears=np.arange(start_year, end_year+1, step=1)
    filler=np.zeros([len(fullyears), len(data_grid.yt_j), len(data_grid.xt_i)])
    P_over_years=xr.DataArray(
        filler, name='P_hum', 
        coords={'Time': fullyears, 'yt_j': data_grid.yt_j, 'xt_i':data_grid.xt_i}, 
        dims=('Time', 'yt_j', 'xt_i'))
    #TEST:
    #P=get_P_map_year(2012, years, data_grid, E_tot, perc, renormalization_factor)        
    #print(P)                  
    
    #The following produces the .nc file and takes roughly 10min
    #import os
    #filename="P_hum_netcdf.nc"
    #if os.path.exists(filename):
    #    os.remove(filename)
    #    print('old file is removed')
    #else:
    #    print("Can not remove the old %s file." % filename)
    import time
    timestr=time.strftime("%h%d_%H-%M")
    txt_files_folder="txt_file_each_year"+timestr+"/"
    os.makedirs(txt_files_folder)
    print("Currently calculating year: ")
    for num,i in enumerate(fullyears):
        if (i%20 == 0): print('%.0f' % i+', ', end='\n')
        if True: print('%.0f' % i+', ', end='')
        
        if i in years:
            P_hum_year=get_P_map_year(i, years, data_grid, E_tot, perc, renormalization_factor)
            P_over_years.sel(Time=i)[:,:]=P_hum_year
        else: #before 1971
            P_hum_year=get_P_map_year(years[0],years, data_grid, E_tot, perc, renormalization_factor)*1.02**(i-1971)
            P_over_years.sel(Time=i)[:,:]=P_hum_year
        txts=open(txt_files_folder+"P_hum_year"+str(i)+".dat", 'w')
        for line in range(0,len(P_hum_year[:,0])):
            for i in range(0,len(P_hum_year[0,:])):
                txts.write(str(P_hum_year[line,i])+"    ")
            txts.write('\n')
        txts.close()
        
    print(P_over_years)
    P_over_years.to_netcdf("P_hum_netcdf_"+timestr+".nc")
    P=xr.open_dataset("P_hum_netcdf_"+timestr+".nc")
    print(P)
    


# Peter Steiglechner, 1.3.
# This program proudces a netcdf file including P_hum 2005 to 2100 for a SSP5 scenario
# The data is taken from the SSP5 and the world regions are somehow adjusted.
#
# Need to implement both versions: 
#   1. in 2010: Jump from WB data to SSP5 data and just take the SSP5 data for future projections
#   2. in 2010: take the WB data and calculate the relative changes from these datasets in the SSP5 scenario.  (then e.g. reforming world ref is overrated!) [so far]
#

from __future__ import print_function
import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import os

regions=["Asia", "Latin America", "Middle East Africa", "OECD", "Reforming", "World"]
r=["as", "lam", "mea", "oecd", "ref", "w"]

def get_Energy():
    data=np.loadtxt("ssp5_marker_world_regions_primE_baseline.txt")
    years=data[0,:]
    E={}
    for n,reg in enumerate(r):
        E[reg]=data[n+1,:]*10**18 /3600/1000 # from EJ to kWh
    return years, E

def get_climber_grid():
    data_grid=xr.open_dataset("snapshots_potsdam2.055623.01.01.dta.nc", decode_times=False)
    return data_grid
        
#### DEFINE NEW AREAS 
# I can mostly use the old areas, but there are a few adjustments...
# E.g. Japan to oecd...        
from create_historic_Phum_files import get_lat_lon_agg
from create_historic_Phum_files import area_lon

def get_lat_lon_agg_ssp(ssp_reg, climber_grid):
    '''
    # Asign each world region to a square (and some additional) cells of CLIMBER.
    # Here not for WorldBank regions but for SSP regions....
    # This builds on get_lat_lon_agg in create_historic_Phum_files
    '''
    array=np.zeros_like(climber_grid)
    minus=[]
    superregion=False
    
    if ssp_reg=='as':
        superregion=True
        array1=get_lat_lon_agg_ssp('eap-aus-jap', climber_grid)  #New region: EAP - Japan and - Australia; defined below
        array2=get_lat_lon_agg('sa', climber_grid) #Old region: South Asia
        array=array1+array2
    if ssp_reg=='eap-aus-jap': # 'Asia':
        lat_low=10; lat_high=17      #same as in PMIP
        lon_low=4; lon_high=6        #same as in PMIP
        additional=[[17,3]]          #same as in PMIP
        minus= [[17,5],[16,5],[17,6], [16,6]]     ###JAPAN

    if ssp_reg=='lam':#'Latin America':
        superregion=True
        array=get_lat_lon_agg('lac', climber_grid)    #Old region: Latin Am and Caribbean

    if ssp_reg=='mea':#'Middle East Africa
        superregion=True
        array1=get_lat_lon_agg('mena', climber_grid)     #Old region: Middle East and north africa
        array2=get_lat_lon_agg('ssa', climber_grid)      #Old region: Subsaharan Africa
        array=array1+array2    
        
    if ssp_reg=='oecd':
        superregion=True
        array1=get_lat_lon_agg('eca', climber_grid)
        array2=get_lat_lon_agg('na', climber_grid)
        array3=get_lat_lon_agg_ssp('jap', climber_grid)    #New region: Japan; defined below
        array4=get_lat_lon_agg('aus', climber_grid)
        array_ref=get_lat_lon_agg_ssp('ref', climber_grid)   #New region: Reforming World (Russia); defined below
        array=array1+array2+array3+array4 - array_ref        # substract the new region: ref

    if ssp_reg=='jap':
        lat_low=16; lat_high=17
        lon_low=5; lon_high=6
        additional=[]
        minus= []    
        
    if ssp_reg=='ref':
        lat_low=17; lat_high=20; lon_low=2; lon_high=2
        additional=[[19,1],[18,1]]; minus=[]
        
    if superregion:
        return array
    else: #still defining new single regions like japan, ref. 
        for lat in range(lat_low, lat_high+1):
            for lon in range(lon_low, lon_high+1):
                array[lat,lon]=1
        for p in additional:
            array[p[0], p[1]]=1
        for p in minus:
            array[p[0], p[1]]=0
    return array
    
    
def get_g(years, full_years, E):
    #from energys in 2005, 2010, ... in ssp, calculate the growth rate (constant) in this period.
    growth_rate={}
    for reg in r:
        growth_rate[reg]=np.zeros([len(years[:-1])])
        for n in range(0,len(years[:-1])):
            growth_rate[reg][n]=(E[reg][n+1]/E[reg][n])**(1./(years[n+1]-years[n]))
            
    # expand the growth rate (still constant throughout the decade) to full_years
    growth_rate_yearly={}
    for reg in r:
        growth_rate_yearly[reg]=np.empty([len(full_years[:-1])])
        for i in range(1,len(years)):
                for k in range(int(years[i-1]), int(years[i])):
                    growth_rate_yearly[reg][k-full_years[0]]=growth_rate[reg][i-1]
    return growth_rate_yearly
    
def get_g_year(year, growth_rate_yearly): #careful when year =2100 no growth rate.
    year_ind=np.where(full_years==year)[0][0]
    g=1.+0*grid.frlnd[0]#np.ones_like(grid.frlnd[0])  #g=1 where there is ocean (and therefore no P)
    for reg in r[:-1]:
        current_area=get_lat_lon_agg_ssp(reg, grid.frlnd[0])
        g+=current_area * (-1+growth_rate_yearly[reg][year_ind])
    return g
    
    
#### FOR OPT 2 ONLY ####
def expand_E_tot_to_fullyears(years, full_years, E_tot):
    # Since E_tot is only defined at 2005, 2010, ... we expand it to every year via interpolation!
    E_tot_yearly={}
    for reg in r:
        E_tot_yearly[reg] = np.interp(full_years, years, E_tot[reg])
        #~ E_tot_yearly[reg]=np.empty([len(full_years[:-1])])
        #~ for i in range(1,len(years)):
                #~ for k in range(int(years[i-1]), int(years[i])):
                    #~ E_tot_yearly[reg][k-2005]=E_tot[reg][i-1]
    return E_tot_yearly
    
    
def get_P_map_year_directSSP5(year, years, data_grid, E_yearly):
    year_ind=np.where(years==year)[0][0]
    P_map=np.zeros_like(data_grid.frlnd[0])
    for num, reg in enumerate(r[:-1]):
        P_yearly_reg=E_yearly[reg][year_ind] *1000/(365*24)   #to get from kWh to W/m^2 
        region_array=get_lat_lon_agg_ssp(reg, data_grid.frlnd[0])
        current_area=region_array*(data_grid.frlnd[0] - data_grid.frlnd[0] * data_grid.frglc[0])
        area_size=area_lon(current_area)*6371000**2
        #P_map += E_tot[reg] * current_area / area_size
        P_map += P_yearly_reg * np.array(region_array) / area_size.values
    return P_map

    
    
if __name__=='__main__':
    timestr=time.strftime("%h%d_%H-%M")
    opt=2
    txt_files_folder="txt_file_each_year"+timestr+"opt"+str(opt)+"/"
    os.makedirs(txt_files_folder)
    
    if opt==1:
        print("Chosen option: start with WB Data in 2010 (NOT 2005) and calculate the changes")
    elif opt==2:
        print("Chosen option: start directly with SSP5 data in 2005. Later i will take all files from 2010 onwards.")
    else:
        print("Choose option!!!")
        exit
    
    start_time=time.time()
    years, E = get_Energy()
    grid=get_climber_grid()
    
    #we will use yearly data not 2005, 2010, 2020, ...
    full_years=np.squeeze(range(int(years[0]), int(years[-1])+1))
    #full_years=full_years[5:]
    print(full_years)
    
    growth_rate_yearly=get_g(years, full_years, E)
    if opt==2:
        E_yearly=expand_E_tot_to_fullyears(years, full_years, E)
    P_historic=xr.open_dataset(climber_path+"data_preparation/PMIP/P_hum/P_hum_netcdf_Jul13_11-20.nc")
    P_10=P_historic.sel(Time=2010).P_hum
    
    filler=np.zeros([len(full_years), len(grid.yt_j), len(grid.xt_i)])
    P_future=xr.DataArray(
        filler, name='P_hum_ssp5', coords={'Time':full_years, 'yt_j':grid.yt_j, 'xt_i':grid.xt_i}, 
        dims=("Time", "yt_j", "xt_i"))
    P_future.sel(Time=2010)[:,:]=P_10
    
    print("Currently calculating year: ")
    for y in full_years[1:]:
        print(str(y)+", ", end='')
        if opt==1:
            #OPTIPON1, i.e. we take WB data and let this grow in each cell according to SSP5
            g_y=get_g_year(y-1, growth_rate_yearly)
            P_hum_year=P_future.sel(Time=y-1)*g_y
            P_future.sel(Time=y)[:,:]=P_hum_year  #P_future.sel(Time=y-1)*g_y
        elif opt ==2:
            #OPTIPON2, i.e. data directly from SSP5
            P_hum_year=get_P_map_year_directSSP5(y, full_years, grid, E_yearly)
            P_future.sel(Time=y)[:,:] = P_hum_year  #get_P_map_year_directSSP5(y, full_years, grid, E_yearly)
        else:
            print("Define opt: 1 == start with WB-data and grow, 2== jump to SSP5")    
    
        txts=open(txt_files_folder+"P_hum_year"+str(y)+".dat", 'w')
        for line in range(0,len(P_hum_year[:,0])):
            for i in range(0,len(P_hum_year[0,:])):
                # CHANGE 13.7. 
                txts.write(str(P_hum_year[line,i])+"    ")
            txts.write('\n')
        txts.close()

    P_future.to_netcdf("P_hum_netcdf_"+timestr+".nc")	
    #P_test=xr.open_dataset("P_hum_netcdf_"+timestr+".nc")
    elapsed_time=time.time()-start_time
    print("Elapsed time = ", elapsed_time, " s")

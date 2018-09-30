import numpy as np

#IMPLEMETNATION
C=4e7	# [J/m^2/K] heat capacity

# parameters
dt=1./4*3600. *1.
nr_years=15
N=int(nr_years*3600*24*365 /dt)
C=4e7

a0=0.3
ai=0.6

#### grid
Nlats=360
lats_star=np.linspace(-90, 90, num=Nlats+1)
x_star=np.sin(lats_star*np.pi/180)
lats=(lats_star[:-1]+ lats_star[1:])/2
x=np.sin(lats*np.pi/180.)
dx=np.diff(x)
dx_star=np.diff(x_star)

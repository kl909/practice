import geopandas as gpd 
import numpy as np
from datetime import timedelta as delta
import xarray as xr
import shapefile
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import PercentFormatter
import pandas as pd


#load dataset
data_xarray = xr.open_zarr('Trajectory_swain_100.zarr')
#data_xarray = xr.open_zarr('Trajectory_test.zarr')

#grab start point coordinates
start_lat = data_xarray["lat"][:,0].values #go to lat and pull all particles, first value
start_lon = data_xarray["lon"][:,0].values                                                                 

#grab end point coordinates
end_lat = data_xarray["lat"][:,-1].values #go to lat and pull all particles, last value                   
end_lon = data_xarray["lon"][:,-1].values                                                                  

#combine start coordinates 
start_points = np.stack([start_lon, start_lat], axis=1) # had to change lat and lon around
end_points = np.stack([end_lon, end_lat], axis=1) # had to change lat and lon around
#print(end_points)

#save startpoints as shapefile
w = shapefile.Writer('start_points')
w.field('name', 'C')
i = 0
for p in start_points:
    w.point(p[0],p[1])
    w.record('point'+str(i))
    i = i+1

w.close()

#save endpoints as shapefile
w = shapefile.Writer('end_points')
w.field('name', 'C')
i = 0
for p in end_points:
    w.point(p[0],p[1])
    w.record('point'+str(i))
    i = i+1

w.close()

################ import in the new shape files

#import end points shapefiles
endf = "end_points.shp"
endpoints = gpd.read_file(endf)

#import start points shp file
startf = "start_points.shp"
startpoints = gpd.read_file(startf)

#import reef polygon shp file
sf = "../GBR_reef_names/Swain_reef/swain_reef.shp"
#sf = "gis_swain_reef/swain_sample/swain_reef_sample.shp"
polys = gpd.read_file(sf)
#print(polys['GBR_NAME'])

################ connectivity matrix
# loop to find which particles start in which reef
start = {}
for p in range(0,len(polys)): # loop over polys
    i = 0
    particle_list = [] # this is the array 
    for pp in range(0,len(startpoints)): # finding which particles start in that polygon
        if (polys['geometry'][p].contains(startpoints['geometry'][pp])):
             particle_list.append(i)
        i = i+1
    start[polys['LOC_NAME_S'][p]] = particle_list #store list in dictionary
#print(start)
#print(endpoints)

# connectivity matrix
matrix = np.zeros((len(polys),len(polys)))
for p in range(0,len(polys)): # loop through polys
    particle_list = start[polys["LOC_NAME_S"][p]] #pulls the particle numbers that started in the poly
    for pp in range(0,len(polys)): #loop over all of the polys again
        for particle in particle_list: #for the particles which started in our "first" polygon
            if (polys['geometry'][pp].contains(endpoints['geometry'][particle])): # which ones finish in our "second" poly
                matrix[p][pp] = matrix[p][pp] + 1 # and add them up
row_names = polys['LOC_NAME_S']
column_names = polys['LOC_NAME_S']
con_matrix = pd.DataFrame(matrix, index=row_names, columns=column_names) #convert to pandas
con_matrix.to_csv('con_matrix.csv') # save to csv

print(con_matrix)

  


df = pd.read_csv("con_matrix.csv")
df = df.drop(columns=['LOC_NAME_S'])

# to make percentage
pc = df.sum() 
pc = (df/pc)*100 
pc=pc.fillna(0.0)

f = plt.figure(figsize=(19, 15))
plt.matshow(pc, fignum=f.number)
cb = plt.colorbar()
cb.ax.set_ylabel('Percentage', size=14)
plt.title('Swain Area', fontsize=18);
plt.savefig('Connectivity_Matrix', bbox_inches='tight')
plt.show()

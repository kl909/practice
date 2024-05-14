import geopandas as gpd 
import numpy as np
from datetime import timedelta as delta
import xarray as xr
import shapefile
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt
import pandas as pd


# section 1 requires firedrake

# section 1: get the coordinates of the trajectories
data_xarray = xr.open_zarr('Trajectory.zarr',decode_times=False)
#make into dataframe
df = data_xarray[['time','lat','lon']].to_dataframe() # don't think I need this with .zarr'
#print(df)

#num of particles and number of steps
numpart = df["time"].last_valid_index()[0]
numstep = df["time"].last_valid_index()[1]

#start point
start_points = []
# make loop to go through the trajectories
for x in range(0, numpart):
   start_points.append([df['lon'][x][0], df['lat'][x][0]])
#print(start_points)

#save list into shapefile
w = shapefile.Writer('start_points')
w.field('name', 'C')
i = 0
for p in start_points:
    w.point(p[0],p[1])
    w.record('point'+str(i))
    i = i+1

w.close()

#section 1 part 2
# extract end points
# make empty list for final coordinates
end_points = []
# make loop to go through the 10 trajectories to give final coordinates
for x in range(0, numpart):
   end_points.append([df['lon'][x][numstep], df['lat'][x][numstep]]) # needed extra set of [] to put lat and lon into list    
#print(end_points)

#save list into shapefile
w = shapefile.Writer('final_points')
w.field('name', 'C')
i = 0
for p in end_points:
    w.point(p[0],p[1])
    w.record('point'+str(i))
    i = i+1

w.close()

#import points shapefiles
endf = "final_points.shp"
endpoints = gpd.read_file(endf)
#print(endpoints)

startf = "start_points.shp"
startpoints = gpd.read_file(startf)
#print(startpoints)
# end of extracting and importing trajectory data


# import reef polygon shapefile
sf = "swain_reef_sample_data.shp"
polys = gpd.read_file(sf)
#print(polys['Reef_name'])



# find starting points of particles --> this needs to be made efficient
start = {}
for p in range(0,len(polys)): # loop over polys
    i = 0
    particle_list = [] # this is the array 
    for pp in range(0,len(startpoints)): # finding which particles start in that polygon
        if (polys['geometry'][p].contains(startpoints['geometry'][pp])):
             particle_list.append(i)
        i = i+1
    start[polys['Reef_name'][p]] = particle_list #store list in dictionary
#print(start)

# matrix of which particles start and end in which polygons
matrix = np.zeros((len(polys),len(polys)))
for p in range(0,len(polys)): # loop through polys
    particle_list = start[polys["Reef_name"][p]] #pulls the particle numbers that started in the poly
    for pp in range(0,len(polys)): #loop over all of the polys again
        for particle in particle_list: #for the particles which started in our "first" polygon
            if (polys['geometry'][pp].contains(endpoints['geometry'][particle])): # which ones finish in our "second" poly
                matrix[p][pp] = matrix[p][pp] + 1 # and add them up
row_names = polys['Reef_name']
column_names = polys['Reef_name']
con_matrix = pd.DataFrame(matrix, index=row_names, columns=column_names) #convert to pandas
print(con_matrix)
#save the matrix
#np.savetxt("larvae_matrix.csv", con_matrix, delimiter=",")
con_matrix.to_csv('con_matrix.csv')

#visualise matrix in matplotlib
plt.matshow(con_matrix)
plt.show()

#data.describe()

import xarray as xr
import geopandas as gpd
import numpy as np
import shapely.wkt
from shapely import geometry
from pathlib import Path
import numpy as np
import shapefile
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon, MultiLineString
import fiona
import pandas as pd
from shapely.geometry import mapping
from pyproj import CRS


# load trajectory zarr data
ds = xr.open_zarr("Trajectory_swain.zarr")

df = ds.to_dataframe()
df['geometry'] = df.apply(lambda x: Point((float(x.lon), float(x.lat))), axis=1)
df = gpd.GeoDataFrame(df, geometry='geometry')
df = df.drop(columns=['time']) # otherwise it won't save
df = df.drop(columns=['z'])
df = df.reset_index() # obs and trajectory were indexes - make them columns
#save to shapefile
df = df.sort_values(by=['trajectory','obs'])

all_trajs = []
for t in range(0,np.max(df["trajectory"])): 
    current_traj = df[df["trajectory"] == t]
    all_trajs.append(LineString(list(current_traj["geometry"])))
    
# Define a feature geometry with one attribute
schema = {
    'geometry': 'LineString',
    'properties': {'id': 'int'},
}

# Write a new Shapefile
with fiona.open('testing123.shp', 'w', crs=CRS.from_epsg(32756).to_proj4(), driver='ESRI Shapefile', schema=schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    count = 0
    for line in all_trajs:
        c.write({
            'geometry': mapping(line), 
            'properties': {'id': count},
           })
        count +=1



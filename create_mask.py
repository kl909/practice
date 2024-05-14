import qmesh3 as qmesh
import sys

input_shape = sys.argv[1]
output_shape = sys.argv[2]

#Initialising qgis API
qmesh.initialise()
#Read in the shapefile describing the domain boundaries and create a polygon of that shape
simuationBoundaries = qmesh.vector.Shapes()
simuationBoundaries.fromFile(input_shape)
loopShapes = qmesh.vector.identifyLoops(simuationBoundaries,
                      isGlobal=False, defaultPhysID=1000,
                      fixOpenLoops=True,
                      extraPointsPerVertex=5)
polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=5e3,
                      smallestMeshedArea=1e10, meshedAreaPhysID = 1)
polygonShapes.writeFile(output_shape)

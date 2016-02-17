import numpy as np

def readVTK(filename, createInterpolant=False, conversionMatrix=[], conversionVector=[], projectionVector=[]):
    """imports standard SOWFA vtk files

    input: file = location of vtk-file
    outputs
    dataType =  OpenFOAM label of measurement (e.g. U, Umean, have not tested for several measurements)
    cellCenters = centers of cells that have been sampled (x,y,z)
    cellData = sampling values (could be vectors (rows))

    to convert to different coordinate frame, you can use conversionMatrix and conversionVector:
    {cellCenter vector new frame} = conversionMatrix * {cellCenter vector VTK frame} +  conversionVector

    to output an interpolator that you can use to sample specific points, use createInterpolant = True

    see example below

    Pieter Gebraad, 2015"""

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    lineCounter = 0

    while lineCounter < len(lines):
        line = lines[lineCounter].strip()

        if line == 'DATASET POLYDATA':
            lineCounter += 1
            line = lines[lineCounter].strip() # read 'POINTS nPoints float'
            line = line.split()
            nPoints = int(line[1])

            pointsXYZ = np.array([ x.split() for x in lines[lineCounter+1:lineCounter+1+nPoints]]).astype(np.float)
            lineCounter = lineCounter+nPoints

        if line[0:8] == 'POLYGONS': # read 'POLYGONS nPolygons ...'
            line = line.split()
            nPolygons = int(line[1])

            polygons = np.array([x.split() for x in lines[lineCounter+1:lineCounter+1+nPolygons]]).astype(np.int)
            polygons = polygons[:, 1:]
            lineCounter = lineCounter+nPolygons

        if line[0:9] == 'CELL_DATA':
            lineCounter += 1
            line = lines[lineCounter].strip()  # read 'FIELD attributes nAttributes'
            line = line.split()
            nAttributes = int(line[2])

            cellData = list()
            dataType = list()
            for att in range(0, nAttributes):
                lineCounter += 1
                line = lines[lineCounter].strip()
                fieldData = line.split()# read 'U 3 nPolygons float'
                if fieldData[3] == 'float' and int(fieldData[2]) == nPolygons:
                    dataType.append(fieldData[0])
                    cd = np.array([ x.split() for x in lines[lineCounter+1:lineCounter+1+nPolygons]]).astype(np.float)
                    if projectionVector != []:
                        cd = np.dot(cd, projectionVector)
                    cellData.append(cd)
                else:
                    print 'readVTK FORMAT ERROR'

        lineCounter += 1

    # cell to point data
    pointsXYZ = pointsXYZ.take(polygons.flatten(1), axis=0)
    cellCenters = np.array([pointsXYZ[:,0].reshape([3,-1]).mean(axis=0),pointsXYZ[:,1].reshape([3,-1]).mean(axis=0),pointsXYZ[:,2].reshape([3,-1]).mean(axis=0)]).transpose()

    if conversionMatrix != [] and conversionVector != []:
        cellCenters = (np.dot(conversionMatrix, cellCenters.transpose()).transpose()) + conversionVector

    if createInterpolant==False:
        if nAttributes == 1:
            dataType = dataType[0]
            cellData = cellData[0]
        return dataType, cellCenters, cellData
    else:
        from scipy.interpolate import griddata
        from copy import copy
        interpolants = list()
        for att in range(0, nAttributes):
            def interpolant(samplePoints):
                sampledData = griddata(cellCenters, cellData[att], samplePoints, method='nearest')
                return sampledData
            interpolants.append(copy(interpolant))
        if nAttributes == 1:
            dataType = dataType[0]
            interpolants = interpolants[0]
        return dataType, interpolants


def averageVTKs(basePath, timeFolders, filename, createInterpolant=True, conversionMatrix=[], conversionVector=[], projectionVector=[]):

    from os import path
    from sys import stdout

    nSamples = len(timeFolders)

    print 'reading and averaging hub-height flow field'
    for dataI in range(nSamples):

        # show progress
        progress = 100.*(dataI+1)/nSamples
        if dataI > 0:
            stdout.write("\b"*16)
        stdout.write("progress: %3s %%" % int(progress))

        dataType, cellCenters, cellData = readVTK(path.join(basePath,timeFolders[dataI],filename), createInterpolant=False, conversionMatrix=conversionMatrix, conversionVector=conversionVector, projectionVector=projectionVector)

        if dataI == 0:
            cellDataMean = cellData/nSamples
        else:
            cellDataMean += cellData/nSamples
    stdout.write("\n")
    if createInterpolant==False:
        return dataType, cellCenters, cellDataMean
    else:
        from scipy.interpolate import griddata
        from copy import copy
        def interpolant(samplePoints):
            sampledData = griddata(cellCenters, cellDataMean, samplePoints, method='nearest')
            return sampledData
        return dataType, copy(interpolant)


if __name__ == '__main__':

    # EXAMPLE
    dataType, cellCenters, cellData = readVTK('U_slice_1.vtk', False)  # replace filename with your own vtk file
    print "type of data: %s" % dataType
    print "cell centers:"
    print cellCenters[0]
    print cellCenters[1]
    print cellCenters[2]
    print "..."
    print cellCenters[-3]
    print cellCenters[-2]
    print cellCenters[-1]
    print "cell data:"
    print cellData[0]
    print cellData[1]
    print cellData[2]
    print "..."
    print cellData[-3]
    print cellData[-2]
    print cellData[-1]

    dataType, interpolant = readVTK('U_slice_1.vtk', True, np.array([[1.0, 0, 0], [0, 1.0, 0]]), np.array([[-1500.0, -1500.0]]))

    import matplotlib.pyplot as plt

    resolution = 600
    #x = np.linspace(-1500, 1500, resolution)
    #y = np.linspace(-1500, 1500, resolution)
    x = np.linspace(0, 3500, resolution*3500.0/3000.0)
    y = np.linspace(0, 3000, resolution)
    xMesh, yMesh = np.meshgrid(x, y)

    velocitiesMesh = interpolant((xMesh.flatten(), yMesh.flatten()))
    absVelocitiesMesh = np.sqrt((velocitiesMesh**2).sum(axis=1)).reshape(resolution,resolution*3500.0/3000.0)

    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.pcolormesh(x, y, absVelocitiesMesh, cmap='coolwarm')
    plt.colorbar(im, orientation='vertical')
    ax1.set_aspect('equal')
    ax1.autoscale(tight=True)
    plt.show()



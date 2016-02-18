from utilities.readVTK import *
import cPickle as pickle

cases = ['neutral', 'stable']
startTimes = [20000, 15000]
casesICOWES3 = ['r01_c1_y30_t-5_R0_I0', 'r03_c3_y30_t-5_R0_I0']

downStreamDiamRange = np.array([1,3,5])
localCopyDir = 'averageProfilesICOWES3'

for caseI, caseICOWES3 in enumerate(casesICOWES3):
    startTime = startTimes[caseI]
    deltaTimeRange = np.arange(705,1000,5)
    timeFolders = startTime + deltaTimeRange
    sdiFolder = '/scratch/pfleming/runs/ICOWES_3/runCases/%s/sliceDataInstant' % caseICOWES3
    fileSlice = 'U_slice_2'
    timeFolders = timeFolders.astype(str)

    dataType, cellCenters, cellData = averageVTKs(sdiFolder, timeFolders, fileSlice+'.vtk', createInterpolant=False)
    d = {'dataType': dataType, 'cellCenters': cellCenters, 'cellData': cellData}
    pickle.dump(d, file(sdiFolder + '/' + fileSlice + '.avg_pickle', 'w'))
    pickle.dump(d, file(localCopyDir + '/' + fileSlice + '.' + caseICOWES3 + '.avg_pickle', 'w'))

    startTime = startTimes[caseI]
    caseICOWES3s = caseICOWES3+'_slices'
    sdiFolder = '/scratch/pfleming/runs/ICOWES_3/runCases/%s/sliceDataInstant' % caseICOWES3s
    deltaTimeRange = np.arange(1000,1365,5)
    timeFolders = startTime + deltaTimeRange
    timeFolders = timeFolders.astype(str)

    for downStreamDiamI, downStreamDiam in enumerate(downStreamDiamRange):
        fileSlice = 'U_slide_%dDdwnstrT1' % downStreamDiam
        dataType, cellCenters, cellData = averageVTKs(sdiFolder, timeFolders,  fileSlice+'.vtk', createInterpolant=False)
        d = {'dataType': dataType, 'cellCenters': cellCenters, 'cellData': cellData}
        pickle.dump(d, file(sdiFolder + '/' + fileSlice + '.avg_pickle', 'w'))
        pickle.dump(d, file(localCopyDir + '/' + fileSlice + '.' + caseICOWES3 + '.avg_pickle', 'w'))
import numpy as np

atomsNumber = np.genfromtxt('scfExtractedData.dat', max_rows = 1)
cellVolume = np.genfromtxt('scfExtractedData.dat', skip_header = 1, max_rows = 1)
latticeParamter = np.genfromtxt('scfExtractedData.dat', skip_header = 2, max_rows = 1)
latticeVectors = np.genfromtxt('scfExtractedData.dat', skip_header = 3, max_rows = 3)
atomsCoordinates = np.genfromtxt('scfExtractedData.dat', skip_header = 6, usecols = (1,2,3))
atomsNames = np.genfromtxt('scfExtractedData.dat', skip_header = 6, usecols = 0, dtype = str)

latticeParamter = latticeParamter*0.52917720859 
latticeParamters = latticeParamter*latticeVectors

print(int(atomsNumber))
print(cellVolume)

for i in range(0, len(latticeParamters)):
    print("{:.6f} {:.6f} {:.6f}".format(latticeParamters[i][0], latticeParamters[i][1], latticeParamters[i][2]))

for i in range(0, len(atomsCoordinates)):
    print(atomsNames[i],"{:.6f} {:.6f} {:.6f}".format(atomsCoordinates[i][0], atomsCoordinates[i][1], atomsCoordinates[i][2]))
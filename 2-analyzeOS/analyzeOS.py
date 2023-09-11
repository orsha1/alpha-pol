import numpy as np
from scipy.spatial import distance
from scipy.constants import value
from collections import Counter

atomsNumber = np.genfromtxt('scfExtractedData.dat', max_rows=1)
cellVolume = np.genfromtxt('scfExtractedData.dat', skip_header=1, max_rows=1)
latticeParameters = np.genfromtxt(
    'scfExtractedData.dat', skip_header=2, max_rows=3)
atomsCoordinates = np.genfromtxt(
    'scfExtractedData.dat', skip_header=5, usecols=(1, 2, 3))
atomsNames = np.genfromtxt('scfExtractedData.dat',
                           skip_header=5, usecols=0, dtype=str)

ASites = ['Na', 'Ba', 'Pb', 'Sr', 'Ca', 'K', 'Sc', 'Bi', 'Sr', 'Li']
ASitesZstar = [1.13, 2.75, 3.92, 2.54, 2.58, 1.14, 'Sc', 'Bi', 2.54, 'Li']
BSites = ['Mg', 'Ni', 'Sb', 'Sn', 'Zn', 'Zr', 'V', 'Ru',
          'Ti', 'Ta', 'W', 'Mo', 'Nb', 'Fe', 'Ga', 'Hf', 'In']
BSitesZstar = [2.0, 'Ni', 'Sb', 'Sn', 'Zn', 6.3, 'V', 4.0,
               6.3, 8.6, 11.0, 11.0, 9.11, 'Fe', 'Ga', 'Hf', 'In']

f = open('aResults.dat', 'w')
fb = open('aShortResults.dat', 'w')
fIG = open('pdfIG.dat', 'w')
fc = open('a.vasp', 'w')

# find A-site and B-site
cellsNumber = int(atomsNumber/5.0)
ASiteCoordinates = np.zeros(shape=(cellsNumber, 3))
BSiteCoordinates = np.zeros(shape=(cellsNumber, 3))
OCoordinates = np.zeros(shape=(int(cellsNumber*3), 3))
a = 0
b = 0
o = 0
ANames = []
BNames = []
AZstars = []
BZstars = []
for i in range(0, len(atomsNames)):
    atom = atomsNames[i]
    if atom in ASites:
        ANames.append(atom)
        AZstars.append(ASitesZstar[ASites.index(atom)])
        for j in range(0, 3):
            ASiteCoordinates[a][j] = atomsCoordinates[i][j]
        a = a + 1
    if atom in BSites:
        BNames.append(atom)
        BZstars.append(BSitesZstar[BSites.index(atom)])
        for j in range(0, 3):
            BSiteCoordinates[b][j] = atomsCoordinates[i][j]
        b = b + 1
    if atom == 'O':
        for j in range(0, 3):
            OCoordinates[o][j] = atomsCoordinates[i][j]
        o = o + 1

# convert to Angstrom
latticeVec = np.zeros(shape=(1, 3))
latticeVec[0][0] = latticeParameters[0][0]
latticeVec[0][1] = latticeParameters[1][1]
latticeVec[0][2] = latticeParameters[2][2]
ASiteCoordinates = ASiteCoordinates*latticeVec
BSiteCoordinates = BSiteCoordinates*latticeVec
OCoordinates = OCoordinates*latticeVec

# Generate extended oxygens network
OExtended = np.zeros(shape=(27*cellsNumber*3, 3))
k = 0
for ix in range(-1, 2):
    for iy in range(-1, 2):
        for iz in range(-1, 2):
            n1 = k*cellsNumber*3
            n2 = (k+1)*cellsNumber*3
            d = latticeVec.copy()
            d[0][0] = d[0][0] * ix
            d[0][1] = d[0][1] * iy
            d[0][2] = d[0][2] * iz
            OExtended[n1:n2] = OCoordinates + d
            k = k+1

# functions for finding neighbors, displacement and alpha


def alphaCalculator(a, b, c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    alpha = 180-np.degrees(np.arccos(cosine_angle))
    if alpha > 180.0:
        alpha = 180-(360.0-alpha)
    return(alpha)


def NeighborsFinder(cation, XO):
    dAll = distance.cdist(cation, XO, 'euclidean')
    i = np.where(dAll <= 3.5)[1]
    d = dAll[0][i]
 #   print(dAll)
 #   print(len(d))
    if len(d) == 6:
        xChain = XO[i].copy()
        yChain = XO[i].copy()
        zChain = XO[i].copy()
        xChain = xChain[(np.abs(xChain[:, 1] - cation[0][1]) <= 0.5)
                        & (np.abs(xChain[:, 2] - cation[0][2]) <= 0.5)]
        yChain = yChain[(np.abs(yChain[:, 0] - cation[0][0]) <= 0.5)
                        & (np.abs(yChain[:, 2] - cation[0][2]) <= 0.5)]
        zChain = zChain[(np.abs(zChain[:, 1] - cation[0][1]) <= 0.5)
                        & (np.abs(zChain[:, 0] - cation[0][0]) <= 0.5)]
        xCM = np.average(xChain, axis=0)
        yCM = np.average(yChain, axis=0)
        zCM = np.average(zChain, axis=0)
        dxCov = (cation[0][0] - xCM[0])
        dyCov = (cation[0][1] - yCM[1])
        dzCov = (cation[0][2] - zCM[2])
        xAlpha = alphaCalculator(xChain[0], cation[0], xChain[1])
        yAlpha = alphaCalculator(yChain[0], cation[0], yChain[1])
        zAlpha = alphaCalculator(zChain[0], cation[0], zChain[1])
        result = "{:.0f} {:.2f} {:.2f} {:.2f} {:.7f} {:.7f} {:.7f}".format(
            len(d), xAlpha, yAlpha, zAlpha, dxCov, dyCov, dzCov)
        shortResults = np.array(
            [[xAlpha, yAlpha, zAlpha, dxCov, dyCov, dzCov]])
    elif len(d) == 12:
        cuboctahedra = XO[i].copy()
        CM = np.average(cuboctahedra, axis=0)
        dx = (cation[0][0] - CM[0])
        dy = (cation[0][1] - CM[1])
        dz = (cation[0][2] - CM[2])
        result = "({:.0f}) {:.7f} {:.7f} {:.7f}".format(len(d), dx, dy, dz)
        shortResults = np.array([[dx, dy, dz]])
    else:
        result = "Number of neighbors ("+str(int(len(d))) + \
            ") found by NeighborsFinder function is not premitted"
    # print(result)
    return(result, XO[i], shortResults)


# use NeighborsFinder to find neighbors
ASiteResults = np.zeros(shape=(cellsNumber, 4))
BSiteResults = np.zeros(shape=(cellsNumber, 7))
displacement = np.zeros(shape=(1, 4))
polarization = np.zeros(shape=(1, 4))
alphaAverage = np.zeros(shape=(1, 3))
# print(len(ASiteCoordinates))
# print(ASiteCoordinates)
# print(len(OExtended))
for i in range(0, len(ASiteCoordinates)):
    res = NeighborsFinder(ASiteCoordinates[i:i+1], OExtended)
    print(ANames[i], AZstars[i], ASiteCoordinates[i:i+1][0],
          (ASiteCoordinates[i:i+1][0]/latticeVec)[0], file=f)
    print(res[0], file=f)
    print(res[1], file=f)
    print(" ", file=f)
    ASiteResults[i][0] = AZstars[i]
    for j in range(0, 3):
        ASiteResults[i][j+1] = res[2][0][j]

for i in range(0, len(BSiteCoordinates)):
    res = NeighborsFinder(BSiteCoordinates[i:i+1], OExtended)
    print(BNames[i], BZstars[i], BSiteCoordinates[i:i+1][0],
          (BSiteCoordinates[i:i+1][0]/latticeVec)[0], file=f)
    print(res[0], file=f)
    print(res[1], file=f)
    print(" ", file=f)
    BSiteResults[i][0] = BZstars[i]
    for j in range(0, 6):
        BSiteResults[i][j+1] = res[2][0][j]

# calculate pols and average alpha
alphaAverage[0][0] = np.average(np.min(BSiteResults[:, 1:4], axis=1))
alphaAverage[0][1] = np.average(np.average(BSiteResults[:, 1:4], axis=1))
alphaAverage[0][2] = np.average(np.max(BSiteResults[:, 1:4], axis=1))

eOmega = value(u'elementary charge')/np.prod((10**(-10))*latticeVec)

BSiteZPx = np.sum(BSiteResults[:, 0]*BSiteResults[:, 4]*(10**(-10)))
ASiteZPx = np.sum(ASiteResults[:, 0]*ASiteResults[:, 1]*(10**(-10)))
BSiteZPy = np.sum(BSiteResults[:, 0]*BSiteResults[:, 5]*(10**(-10)))
ASiteZPy = np.sum(ASiteResults[:, 0]*ASiteResults[:, 2]*(10**(-10)))
BSiteZPz = np.sum(BSiteResults[:, 0]*BSiteResults[:, 6]*(10**(-10)))
ASiteZPz = np.sum(ASiteResults[:, 0]*ASiteResults[:, 3]*(10**(-10)))

polarization[0][0] = eOmega*(BSiteZPx+ASiteZPx)
polarization[0][1] = eOmega*(BSiteZPy+ASiteZPy)
polarization[0][2] = eOmega*(BSiteZPz+ASiteZPz)
polarization[0][3] = np.sqrt(
    polarization[0][0]**2 + polarization[0][1]**2 + polarization[0][2]**2)
print("{} {:.7f} {:.7f} {:.7f} {:.7f} {:.7f} {:.7f} {:.7f}".format("outOS", polarization[0][0], polarization[0][1], polarization[0][2],
                                                                   polarization[0][3],
                                                                   alphaAverage[0][0], alphaAverage[0][1], alphaAverage[0][2]), file=fb)
# generate input for IG
AList = list(set(ANames))
BList = list(set(BNames))
print(cellsNumber, len(AList), len(BList), "1",
      "num_uc num_spec_a num_spec_b num_spec_x", file=fIG)

for i in range(0, len(AList)):
    print(ANames.count(AList[i]), AList[i],
          "    0.9 0.045 4.5 num_ions name bscat therm", file=fIG)

for i in range(0, len(BList)):
    print(BNames.count(BList[i]), BList[i],
          "    0.9 0.045 4.5 num_ions name bscat therm", file=fIG)
print(cellsNumber*3, "O", "    0.9 0.045 4.5 num_ions name ascat therm", file=fIG)

for i in range(0, len(AList+BList)):
    print("1.904 6.0 R N", file=fIG)

print(latticeVec[0][0], latticeVec[0][1], latticeVec[0][2], file=fIG)

for j in range(0, len(AList)):
    for i in range(0, len(atomsNames)):
        if AList[j] == atomsNames[i]:
            print(atomsCoordinates[i][0], atomsCoordinates[i]
                  [1], atomsCoordinates[i][2], file=fIG)

for j in range(0, len(BList)):
    for i in range(0, len(atomsNames)):
        if BList[j] == atomsNames[i]:
            print(atomsCoordinates[i][0], atomsCoordinates[i]
                  [1], atomsCoordinates[i][2], file=fIG)

for i in range(0, len(atomsNames)):
    if atomsNames[i] == 'O':
        print(atomsCoordinates[i][0], atomsCoordinates[i]
              [1], atomsCoordinates[i][2], file=fIG)

# generate vasp for vesta
lB = [Counter(BNames)[BList[i]] for i in range(len(BList))]
lA = [Counter(ANames)[AList[i]] for i in range(len(AList))]
print("New structure", file=fc)
print("1.0", file=fc)
print("        {:.7f} {:.7f} {:.7f}".format(
    latticeParameters[0][0], latticeParameters[0][1], latticeParameters[0][2]), file=fc)
print("        {:.7f} {:.7f} {:.7f}".format(
    latticeParameters[1][0], latticeParameters[1][1], latticeParameters[1][2]), file=fc)
print("        {:.7f} {:.7f} {:.7f}".format(
    latticeParameters[2][0], latticeParameters[2][1], latticeParameters[2][2]), file=fc)
print(*BList, *AList, " O", sep=' ', file=fc)
print(*lB, *lA, cellsNumber*3, file=fc)
print("Direct", file=fc)
for j in range(0, len(BList)):
    for i in range(0, len(atomsNames)):
        if atomsNames[i] == BList[j]:
            print("    {:.7f} {:.7f} {:.7f}".format(
                atomsCoordinates[i][0], atomsCoordinates[i][1], atomsCoordinates[i][2]), file=fc)
for j in range(0, len(AList)):
    for i in range(0, len(atomsNames)):
        if atomsNames[i] == AList[j]:
            print("    {:.7f} {:.7f} {:.7f}".format(
                atomsCoordinates[i][0], atomsCoordinates[i][1], atomsCoordinates[i][2]), file=fc)
for i in range(0, len(atomsNames)):
    if atomsNames[i] == "O":
        print("    {:.7f} {:.7f} {:.7f}".format(
            atomsCoordinates[i][0], atomsCoordinates[i][1], atomsCoordinates[i][2]), file=fc)

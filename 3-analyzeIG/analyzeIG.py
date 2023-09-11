import numpy as np
from scipy.constants import value

ASites = ['Na', 'Ba', 'Pb', 'Sr', 'Ca', 'Kk', 'Sc', 'Bi', 'Sr', 'Li']
ASitesZstar = [1.13, 2.75, 3.92, 2.54, 2.58, 1.14, 'Sc', 'Bi', 2.54, 'Li']
BSites = ['Mg', 'Ni', 'Sb', 'Sn', 'Zn', 'Zr', 'Vv', 'Ru', 'Ti', 'Ta', 'Ww', 'Mo', 'Nb', 'Fe', 'Ga', 'Hf', 'In']
BSitesZstar = [2.0, 'Ni', 'Sb', 'Sn', 'Zn', 6.0, 'Vv', 4.0, 7.06, 8.6, 11.0, 11.0, 9.11, 'Fe', 'Ga', 'Hf', 'In']

f = open('bResults.dat', 'w')
fb = open('bShortResults.dat', 'w')
fc = open('bTiltResults.dat', 'w')
fd = open('bTiltShortResults.dat', 'w')


dispCov = np.genfromtxt('outDispCov.dat', usecols = (2,3,4,5))
dispCovNames = np.genfromtxt('outDispCov.dat', usecols = 0, dtype = str)
dispCovNames = [s.replace('-O_complex_data:' , '') for s in dispCovNames]

disp = np.genfromtxt('outDisp.dat', usecols = (2,3,4,5))
dispNames = np.genfromtxt('outDisp.dat', usecols = 0, dtype = str)
dispNames = [s.replace('-O_complex_data:' , '') for s in dispNames]

tilt = np.genfromtxt('outTilt.dat', usecols = (2,3,4))
tiltNames = np.genfromtxt('outTilt.dat', usecols = 0, dtype = str)
tiltNames = [s.replace('-O_complex_data:' , '') for s in dispNames]

BSitesUZstar = dispCov.copy()
ASitesUZstar = disp.copy()

eOmega = value(u'elementary charge')/(len(BSitesUZstar)*64*(10**(-30)))
n = 0
for i in range(0, len(dispCovNames)):
    if dispCovNames[i] in BSites:
        if isinstance(BSitesZstar[BSites.index(dispCovNames[i])],float):
            BSitesUZstar[n] = BSitesZstar[BSites.index(dispCovNames[i])]*dispCov[i]*(10**(-10))
            n = n + 1
        else:
            print("Error, No ZStar for "+dispCovNames[i])
    else:
        print("Error, "+dispCovNames[i]+" is not in BSite")
n = 0
for i in range(0, len(dispNames)):
    if dispNames[i] in ASites:
        if isinstance(ASitesZstar[ASites.index(dispNames[i])],float):
            ASitesUZstar[n] = ASitesZstar[ASites.index(dispNames[i])]*disp[i]*(10**(-10))
            n = n + 1
        else:
            print("Error, No ZStar for "+dispNames[i])
    else:
        print("Error, "+dispNames[i]+" is not in ASite")

polarization = np.zeros(shape = (1,4))
tiltTot = np.zeros(shape = (1,4))

polarization[0][0] = eOmega*(np.sum(BSitesUZstar[:,0]) + np.sum(ASitesUZstar[:,0]))
polarization[0][1] = eOmega*(np.sum(BSitesUZstar[:,1]) + np.sum(ASitesUZstar[:,1]))
polarization[0][2] = eOmega*(np.sum(BSitesUZstar[:,2]) + np.sum(ASitesUZstar[:,2]))
polarization[0][3] = np.sqrt(polarization[0][0]**2 + polarization[0][1]**2 + polarization[0][3]**2) 

tiltTot[:,0] = np.sum(tilt[:,0], axis = 0)
tiltTot[:,1] = np.sum(tilt[:,1], axis = 0)
tiltTot[:,2] = np.sum(tilt[:,2], axis = 0)
tiltTot[:,3] = tiltTot[0][0] + tiltTot[0][1] + tiltTot[0][3]


print("{} {:.7f} {:.7f} {:.7f} {:.7f}".format("outIG", polarization[0][0], polarization[0][1], polarization[0][2], polarization[0][3]), file = fb)
print("{} {:.7f} {:.7f} {:.7f} {:.7f}".format("tilt", tiltTot[0][0], tiltTot[0][1], tiltTot[0][2], tiltTot[0][3]), file = fd)

for i in range(0, len(tiltNames)):
    print("{:.7f} {:.7f} {:.7f}".format(tilt[i][0], tilt[i][1], tilt[i][2]), file = fc)

for i in range(0, len(dispNames)):
    print("{}: ({:.7f}  {:.7f}  {:.7f}  {:.7f}) ({:.7f}  {:.7f}  {:.7f}  {:.7f})".format(dispNames[i], disp[i][0], disp[i][1], disp[i][2], disp[i][3],
                                                                                         ASitesUZstar[i][0]*(10**(10)), ASitesUZstar[i][1]*(10**(10)), ASitesUZstar[i][2]*(10**(10)), ASitesUZstar[i][3]*(10**(10))), file=f)
for i in range(0, len(dispCovNames)):
    print("{}: ({:.7f}  {:.7f}  {:.7f}  {:.7f}) ({:.7f}  {:.7f}  {:.7f}  {:.7f})".format(dispCovNames[i], dispCov[i][0], dispCov[i][1], dispCov[i][2], dispCov[i][3],
                                                                                         BSitesUZstar[i][0]*(10**(10)), BSitesUZstar[i][1]*(10**(10)), BSitesUZstar[i][2]*(10**(10)), BSitesUZstar[i][3]*(10**(10))), file=f)


sed -i 's/K/Kk/g' Input.dat
sed -i 's/V/Vv/g' Input.dat
sed -i 's/W/Ww/g' Input.dat
gfortran newpdf4.f
./a.out
grep DispCov *O* > outDispCov.dat
grep Disp *O* | grep -v DispCov | grep -v DispTot> outDisp.dat
grep Tilt *O* > outTilt.dat
sed -i 's/\*//g' outTilt.dat
sed -i 's/\*//g' outDisp.dat
sed -i 's/\*//g' outDispCov.dat

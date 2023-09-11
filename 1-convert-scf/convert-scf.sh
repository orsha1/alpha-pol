fileName="scf.out" # input file name

atomsNumber=$(grep "number of atoms" $fileName | awk '{print $5}') # number of atoms

latticeVector1=$(grep "a(1)" $fileName | awk '{print $4, $5, $6}') # lattice vectors
latticeVector2=$(grep "a(2)" $fileName | awk '{print $4, $5, $6}')
latticeVector3=$(grep "a(3)" $fileName | awk '{print $4, $5, $6}')

latticeParamter=$(grep "celldm(1)" $fileName | awk '{print $2}') # lattice parameter
cellVolume=$(grep "unit-cell volume" $fileName | awk '{print $4}') # unit-cell volume

echo $atomsNumber > scfExtractedData.dat
echo $cellVolume >> scfExtractedData.dat
echo $latticeParamter >> scfExtractedData.dat
echo $latticeVector1 >> scfExtractedData.dat
echo $latticeVector2 >> scfExtractedData.dat
echo $latticeVector3 >> scfExtractedData.dat
grep -A$atomsNumber "positions (cryst. coord.)" $fileName | awk '{print $2, $7, $8, $9}' | grep -v "n." >> scfExtractedData.dat # atoms coordinates 

python3 convert-scf.py > scfExtractedData2.dat
mv scfExtractedData2.dat scfExtractedData.dat
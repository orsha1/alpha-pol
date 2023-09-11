bash clear.sh
for fname in BZT_40atoms_2x2x2sys_vcrelax_dis100_vdw-df-ob86.opt #001-001-1-scf.out  001-001-2-scf.out  001-001-3-scf.out  001-001-4-scf.out  001-001-5-scf.out  111-001-1-scf.out  111-001-2-scf.out  111-001-3-scf.out  111-001-4-scf.out  111-001-5-scf.out
do
    echo "starting with ${fname}"
    cp ./0-scf/$fname ./1-convert-scf/scf.out
    cd ./1-convert-scf/
    bash convert-scf.sh
    cd ../2-analyzeOS/
    cp ../1-convert-scf/scfExtractedData.dat .
    python3 analyzeOS.py
    cd ../3-analyzeIG/
    #find . -maxdepth 1 ! -name 'newpdf4.f' ! -name 'analyzeIG.sh' ! -name 'analyzeIG.py' -type f -exec rm {} +
    cp ../2-analyzeOS/pdfIG.dat Input.dat
    #gfortran newpdf4.f 
    #./a.out
    bash analyzeIG.sh
    python3 analyzeIG.py
    cd ../4-outputs
    rm -rf ${fname}
    mkdir ${fname}
    cp ../0-scf/${fname} ./${fname}/
    cp ../1-convert-scf/scfExtractedData.dat ./${fname}/${fname}_ExtractedData.dat
    cp ../2-analyzeOS/aResults.dat ./${fname}/${fname}_ResultsOS.dat
    cp ../2-analyzeOS/aShortResults.dat ./${fname}/${fname}_ShortResultsOS.dat
    cp ../2-analyzeOS/a.vasp ./${fname}/${fname}.vasp
    cp ../3-analyzeIG/Input.dat ./${fname}/${fname}_InputIG.dat
    cp ../3-analyzeIG/Input.dat ./${fname}/${fname}_InputIG.dat
    cp ../3-analyzeIG/bResults.dat ./${fname}/${fname}_ResultsIG.dat
    cp ../3-analyzeIG/bShortResults.dat ./${fname}/${fname}_ShortResultsIG.dat
    cp ../3-analyzeIG/bTiltResults.dat ./${fname}/${fname}_ResultsTilt.dat
    cp ../3-analyzeIG/bTiltShortResults.dat ./${fname}/${fname}_ShortResultsTilt.dat
    cd ../
done
bash clear.sh
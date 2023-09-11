rm -rf 1-convert-scf/*out
rm -rf 1-convert-scf/*dat
rm -rf 2-analyzeOS/*dat
cd 3-analyzeIG/
find . -maxdepth 1 ! -name 'newpdf4.f' ! -name 'analyzeIG.sh' ! -name 'analyzeIG.py' -type f -exec rm {} +

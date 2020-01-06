# Requirements: Python 3.7 with numpy and matplotlib

sudo apt install build-essential gfortran g++
git clone https://github.com/estherlin/dna-origami

curl http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz -so mfold-3.6.tar.gz
tar -xvzf mfold-3.6.tar.gz
cp dna-origami/mfold_quik mfold-3.6/scripts/mfold_quik.in

cd mfold-3.6/
./configure
make
sudo make install

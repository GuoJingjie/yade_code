#! /bin/bash

# use whichever exists. ver.6 on ubuntu 18.04, ver.8 on ubuntu 20.04
ls -la /root/OpenFOAM/OpenFOAM-6/etc/bashrc && source  /root/OpenFOAM/OpenFOAM-6/etc/bashrc
ls -la /root/OpenFOAM/OpenFOAM-8/etc/bashrc && source  /root/OpenFOAM/OpenFOAM-8/etc/bashrc

cd ..
rm -rf Yade-OpenFOAM-coupling
git clone https://github.com/dpkn31/Yade-OpenFOAM-coupling.git
cd Yade-OpenFOAM-coupling
git checkout yadeTestPar
./Allclean
./Allwmake

cd ../trunk/examples/openfoam/example_icoFoamYade
echo `pwd`
echo "---------- creating symlink yadeimport.py ----------"
ln -s ../../../install/bin/yade-ci ./yadeimport.py
ls -la ./yadeimport.py
blockMesh
decomposePar
mkdir yadep
mkdir spheres

mpiexec --allow-run-as-root -n 1 python3 scriptYade.py : -n 2 icoFoamYade -parallel
mpiexec --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py : -n 2 icoFoamYade -parallel


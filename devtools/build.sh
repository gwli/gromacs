TOP=`pwd`/..
installDir=$TOP/gromacs-install
cd $TOP
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=$installDir
make
make check
sudo make install
source $installDir/bin/GMXRC

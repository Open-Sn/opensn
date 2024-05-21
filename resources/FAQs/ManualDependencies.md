
## Manual installation of PETSc
```shell
./configure  \
--prefix=$PWD/../petsc-3.17.0-install  \
--download-hypre=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-shared-libraries=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--with-cxx-dialect=C++11  \
--with-64-bit-indices \
CC=$CC CXX=$CXX FC=$FC \
CFLAGS='-fPIC -fopenmp'  \
CXXFLAGS='-fPIC -fopenmp'  \
FFLAGS='-fPIC -fopenmp'  \
FCFLAGS='-fPIC -fopenmp'  \
F90FLAGS='-fPIC -fopenmp'  \
F77FLAGS='-fPIC -fopenmp'  \
COPTFLAGS='-O3 -march=native -mtune=native'  \
CXXOPTFLAGS='-O3 -march=native -mtune=native'  \
FOPTFLAGS='-O3 -march=native -mtune=native'  \
PETSC_DIR=$PWD
```

## Manual installation of VTK
In the `VTK` directory do
```shell
export VTK_INSTALL_DIR=$PWD/vtk-9.1.0-install
```
now navigate to the extracted `vtk-9.1.0` directory.
```shell
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=$VTK_INSTALL_DIR \
-DBUILD_SHARED_LIBS:BOOL=ON \
-DVTK_Group_MPI:BOOL=ON \
-DVTK_GROUP_ENABLE_Qt=NO \
-DVTK_GROUP_ENABLE_Rendering=NO \
-DVTK_GROUP_ENABLE_Imaging=NO \
-DVTK_GROUP_ENABLE_StandAlone=WANT \
-DVTK_GROUP_ENABLE_Web=NO \
-DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
-DVTK_BUILD_TESTING:BOOL=OFF \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS=-std=c++11 \
../
```
Then make
```shell
make -j8
```
Then install
```shell
make install
```


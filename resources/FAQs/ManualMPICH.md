# Compiling MPICH from source on unix machines

## Step 1
Download the `tar.gz` of the version you want from the mirror site you want.
The untar it, for example
```bash
tar -zxf mpich-4.0.2.tar.gz
cd mpich-4.0.2.tar.gz
```

## Step 2
Make a build directory
```bash
mkdir gcc-build
cd gcc-build
```

## Step 3
Configure the build. Use `--prefix=<where-u-want-it-installed>` to set the
eventual output of `make install` (which comes later). This prefix is where the
`bin`, `include`, `lib`, etc. folder will be. In my example below I want it to
be somewhere in my home directory.

```bash
../configure --prefix=$PWD/../gcc-install \
--enable-shared \
--enable-sharedlibs=gcc \
--enable-fast=O2 \
--enable-debuginfo \
--enable-totalview \
--enable-two-level-namespace \
CC=gcc \
CXX=g++ \
FC=gfortran \
F77=gfortran \
F90='' \
CFLAGS='' \
CXXFLAGS='' \
FFLAGS='-fallow-argument-mismatch' \
FCFLAGS='-fallow-argument-mismatch' \
F90FLAGS='' \
F77FLAGS=''
```

## Step 4
Make it. 
```bash
make -j4
```

## Step 6
Install it.
```bash
make install
```
If you are an admin you probably would want to wrap this in sudo.
# Compiling gcc from source on linux machines

## Step 1
Download the `tar.gz` of the version you want from the mirror site you want.
The untar it, for example
```bash
tar -zxf gcc-10.2.0.tar.gz
cd gcc-10.2.0
```

## Step 2
Download the prerequisites
```bash
./contrib/download_prerequisites
```

## Step 3
Make a build directory
```bash
mkdir gcc-build
cd gcc-build
```

## Step 4
Configure the build. Use `--prefix=<where-u-want-it-installed>` to set the 
eventual output of `make install` (which comes later). This prefix is where the
`bin`, `include`, `lib`, etc. folder will be. In my example below I want it to
be somewhere in my home directory.

```bash
../configure --prefix=~/software/GCC/gcc-10.2.0 \
--disable-multilib \
--enable-languages=c,c++,fortran,jit \
--enable-checking=release \
--enable-host-shared \
--with-pic
```

When on MacOS: if the install complains about something related to 
`/usr/install` then find the latest MacOS SDK, e.g., 
`/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/` and pass this as 
an additional argument:
```bash
../configure --prefix=$PWD/../gcc-install \
--disable-multilib \
--enable-languages=c,c++,fortran,jit \
--enable-checking=release \
--enable-host-shared \
--with-pic \
--with-sysroot=/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/
```

## Step 5
Make it. Grab a cup of coffee this takes a long time. I used `j16` and it took 
14 minutes.
```bash
make -j4
```

## Step 6
Install it.
```bash
make install
```
If you are an admin you probably would want to wrap this in sudo.

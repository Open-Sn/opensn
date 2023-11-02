# Easy install on MacOS Machines

The following instructions were tested on MacOS Monterey (= version 12).

### Step 0

Install ```Homebrew``` (Mac ports could be an alternative).

### Step 1 - Installing GCC, GFortran and the basic environment

Then, you may need to install, using brew, the following:
- make
- cmake
- gcc
- wget

Each time, use: 
```bash 
brew install *package_name*
```

I created a ```/Users/USERNAME/local/bin/``` directory where I symbolically linked 
```gcc```, ```g++```, ```gfortran```
```bash
cd /Users/USERNAME/local/bin/
ln -s /Users/USERNAME/local/homebrew/Cellar/gcc/12.2.0/bin/gcc-12 gcc
ln -s /Users/USERNAME/local/homebrew/Cellar/gcc/12.2.0/bin/gfortran-12 gfortran
ln -s /Users/USERNAME/local/homebrew/Cellar/gcc/12.2.0/bin/g++-12 g++
```

In your ```.bashrc``` file (or ```.bash_profile```), add these lines_
```bash
export PATH=/Users/jean.ragusa/local/bin/:$PATH

export CC=/Users/USERNAME/local/local/bin/gcc
export CXX=/Users/USERNAME/local/local/bin/g++
export FC=/Users/USERNAME/local/local/bin/gfortran
```
Source your ```.bashrc``` to take into accounts those updates:
```bash
source ~/.bashrc
```
Make sure the following environment variables have been properly set
```bash
$CC --version
$CXX --version
$FC --version
```
### Step 2 - MPICH

Download ```MPICH``` from the web locally to your drive. Untar the ball. 

Configure and install as follows:
```bash
./configure --prefix=$PWD/build FFLAGS=-fallow-argument-mismatch FCFLAGS=-fallow-argument-mismatch
make
make install
```

Check that everything went fine:
```bash
mpicc --version
```
Now, update your ```.bashrc``` file
```bash
export PATH=/Users/USERNAME/local/MPICH/mpich-4.0.3/build/bin:$PATH
 
export CC=/Users/USERNAME/local/MPICH/mpich-4.0.3/build/bin/mpicc
export CXX=/Users/USERNAME/local/MPICH/mpich-4.0.3/build/bin/mpicxx
export FC=/Users/USERNAME/local/MPICH/mpich-4.0.3/build/bin/mpifort
```
which means that, from now on, the C, C++, and Fortran compilers are the ones from MPICH.

Make sure to ```source ~/.bashrc``` (or ```source ~/.bash_profile```) to take into accounts those updates.

### Step 3 - Clone Chi-Tech

**Note:** From now on, the rest of the instructions are identical to the ones from the 
[Easy Linux instructions](./Install_ubuntu_easy.md)

**Important:**  If you want to contribute to **Chi-Tech**, it is strongly recommended to first fork the **Chi-Tech** repository into your own Git account and then to clone your fork. 

Clone the **Chi-Tech** repository.  Go the folder where you want to keep Chi-Tech relevant stuff:
```bash
    $ git clone https://github.com/chi-tech/chi-tech
```
or
```bash
    $ git clone https://github.com/YOUR-NAME/chi-tech
```

**Important:** We recommend building all the dependencies into a separate folder. For instance_,
```bash
    $ mkdir dependencies
```

Before building the dependencies, you need to export a few variables for the PETSc install:
```bash
    $ export CC=mpicc
    $ export CXX=mpicxx
    $ export FC=mpifort
```

Go to the chi-tech folder you have just cloned and type:
```bash
    $ cd chi-tech
    $ python3 resources/configure_dependencies.py -d ../dependencies
```
The configure script will attempt to download and install all the necessary 
dependencies **and may take a long time**

### Step 4 - Configure environment

The next step in this process is to setup the environment variables for compiling
Chi-Tech.

```bash
    $source ../dependencies/configure_deproots.sh
```
**Note:** You can replace ```$source ``` in the above with ```$. ```

### Step 5 - Build Chi-Tech

To compile Chi-Tech now just execute:
```bash
    $ ./configure.sh
```
The configure script will generate the CMake build scripts.

In the main directory (i.e. *chi-tech/*), execute:
```bash
    $ make -j4
```

You can also use -j8 even if you don't have 8 processors, the make command
will use threading where possible.

Whenever recompilation is needed or the configuration has changed,
you can run
```bash
    $ ./configure.sh clean
    $ ./configure.sh
    $ make -j4
```

### Step 6 - Run regression tests

To check if the code compiled correctly, execute the test scripts:

```bash
    $ test/run_tests -d test/ -j8
```

### Step 8 - Chi-Tech documentation

You can either access the documentation online [here](https://chi-tech.github.io), or generate it locally.

To generate the documentation from your local working copy, first make sure
Doxygen and LaTeX are installed:

```bash
sudo apt-get install doxygen texlive
```

The documentation is contained in the *doc* folder and can be generated
using a script provided in that folder:

```bash
cd doc/
./YReGenerateDocumentation.sh
```

Once finished, you can view the generated documentation by opening

```bash
doc/HTMLdocs/html/index.html
```

in your favorite browser.
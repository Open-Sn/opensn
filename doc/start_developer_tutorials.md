# Developer Tutorials

## Lua for developers

When developing a new application using Chi-Tech, you need to register some objects to be 
accessible from the lua console. We consider two situations:
1. you are using Chi-Tech as a library
2. you are developing directly in ```modules```


## Tutorials

The developer tutorials are accessible on the documentation, by clicking on the left column under
[Developer's manual](https://chi-tech.github.io/d7/db6/_programmer_manual.html).

## A fully-detailed example:
This is taken from my [chi-ip](https://github.com/ragusa/chi-ip), where I added
IP for DFEM diffusion in Chi-Tech. 

### Step-0:

1. Create empty repo on github
2. git clone (with token, if you want)
3. Open new project in CLion and point to repo that was jsut cloned

### Step-1:
Follow the instructions in the developer's manual [Using Chi-Tech as a library](https://chi-tech.github.io/d8/d99/_dev_man_using_lib.html).

Specifically,
- create the CMakeLists.txt file and populate
- create the main .cc file (```testIP.cc```, in my example) and populate
- Go to CLion's Cmake setting and add the ENV variable for where the amin chi-tech sources are 
(for me, ```CHI_TECH_DIR=/Users/jean.ragusa/repo/chi-tech```)
    - however, this was not sufficient when I compile from a terminal outside of CLion. I had to export that env 
  variable also in that terminal
- check in the CMake terminal that everything went fine (that Downstream.cmake was found, that MPI was found, ...)

### Step-2:

Start adding source files.
- I created a new folder (```DFEMDIffusionSolver```), added source code there.
- I had to update the ```CMakeLists.txt``` file at the upper level by adding ```add_subdirectory("${PROJECT_SOURCE_DIR}/DFEMDiffusionSolver/")```
- I created a new ```CMakeLists.txt``` file in the new folder and included the following lines_:
   ```bash
   file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
   set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE) 
   ```

### Step-3:

- I developed new code in the new folder ```DFEMDIffusionSolver```.
- I add a subfolder there, called ```lua``` where the lua wrappers are found. The wrappers allow to call the functions 
from the lua console (i.e., the input file)
    - ```create_solver.cc```, the solver itself
    - ```setBCproperty.cc```, the BC conditions
    - ```disc_ord_steady_state_lua_utils.h```, the header
- the main code file, ```testIP.cc```, needs to register the lua functions

Later, when the code gets folded into the main chi-tech, a few changes will need to take place. We describe those later.

### Step-4:

I wanted to create a simulation test that uses the new solver, but also contains a new function (the L2 error norm) 
that I do not want to be in the main solver.

- I created a new folder structure that mimics the Chi-Tech simulation test folder structure: ```framework\LuaTest\ ```
- There, I added ```sim_IP_MMS_L2error.cc```, that basically calls my original solver, but also contains the L2 error 
computation. This file is actually the lua wrapper, same functionality as above.
- I also added a header file, ```unit_tests.h```.
- A ```CmakeLists.txt``` file needs to be added there
  ```bash
  file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
  set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
   ```
- In the main code test, ```testIP.cc```, the new lua function needs to be registered.
- In the main ```CmakeLists.txt``` file, the new subfolder needs to be added ```add_subdirectory("${PROJECT_SOURCE_DIR}/framework/LuaTest/") ```

### Step-5:

When bringing the code into the main chi-tech repository,
the following changes need to happen:
- The folder ```DFEMDIffusionSolver``` gets placed in
  ```modules``` and the ```disc_ord_steady_state_lua_utils.h``` in that subfolder
  receives an update: we add ```RegisterLuaEntities()``` inside a namespace
- We add
    - this line ```#include "DFEMDiffusion/lua/ip_lua_utils.h" ```
    - and this line ```dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L)```

  in ```modules\lua\chi_modules_lua.cc```
- Add the proper subfolder in the ```CMakeLists.txt``` file
  located in ```modules``` (i.e., add this line ```add_subdirectory("DFEMDiffusion")```)
- Update ```lua_test.cc``` to register the sim test


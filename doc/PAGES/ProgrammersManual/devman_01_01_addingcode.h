/**\page DevManAddingCode Adding code to the system

The underlying build system used by ChiTech is cmake and its entry point
is the *CMakeLists.txt* file located in the root directory.

\section devman0_sec0 So you want to add your own code

\subsection devman0_sec0_1 Step 1 - Create your directory

Create a directory where you want to add your code. For example, let us
 suppose you want to add a new module. Add a directory to the
 CHI_TECH/CHI_MODULES folder:

\image html devman_newdir.gif "Creating a new directory" width=600px

\subsection devman0_sec0_2 Step 2 - If there is no CMakeLists.txt, create one

In the folder you just created, make a file called "CMakeLists.txt" and
add the following code to it

\code
file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
\endcode

This code will get executed by cmake and will recursively find all files
 with a .cc extension and compile it.

\subsection devman0_sec0_3 Step 3 - Create your source code

Create appropriate headers and source code.

\subsection devman0_sec0_4 Step 4 - Add the folder to the master CMakeLists.txt

The final step is to add the folder you created to the master *CMakeLists.txt*
 document contained in the root folder. If the new folder is already recursed by
 a higher level *CMakeLists.txt* file then this step is not needed.

Look for the line with comment "Define source directories" and your folder
to the "add_subdirectory" logic:

\code
#================================================ Define source directories
set(SOURCES "${CHI_TECH_DIR}/chi_tech_main.cc")
add_subdirectory("${CHI_TECH_DIR}/CHI_CONSOLE")
add_subdirectory("${CHI_TECH_DIR}/CHI_LIB")
add_subdirectory("${CHI_TECH_DIR}/CHI_LUA")
add_subdirectory("${CHI_TECH_DIR}/CHI_MATH")
add_subdirectory("${CHI_TECH_DIR}/CHI_PHYSICS")
add_subdirectory("${CHI_TECH_DIR}/CHI_GRAPH")

add_subdirectory("${CHI_TECH_DIR}/CHI_TIMER")
add_subdirectory("${CHI_TECH_DIR}/CHI_TOOLS")
add_subdirectory("${CHI_TECH_DIR}/CHI_MESH")
add_subdirectory("${CHI_TECH_DIR}/CHI_MPI")
add_subdirectory("${CHI_TECH_DIR}/chi::log")

add_subdirectory("${CHI_TECH_MOD}/CHI_MONTECARLON")
add_subdirectory("${CHI_TECH_MOD}/CHI_DIFFUSION")
add_subdirectory("${CHI_TECH_MOD}/CHI_NPTRANSPORT")

add_subdirectory("${CHI_TECH_MOD}/MyTestModule")
\endcode

\subsection devman0_sec0_5 Step 5 - Include headers and use the code

That's it! The cmake system needs to be run again by executing the
configure script.

\code
./configure.sh
\endcode

And your code will be linked in. Now just make as usual.

\code
make -j5
\endcode

\date Aug 19, 2019

 */
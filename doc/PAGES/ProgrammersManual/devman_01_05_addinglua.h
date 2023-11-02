/**\page DevManAddingLua Adding Lua-input to the system

Lua is the main scripting system used by ChiTech.

\section devman03_sec0 Adding Lua-input

## Step 1: Create a source file
First create a `.cc` file in the location where you want the lua-input to
reside. Example: `CHI_TECH/LuaTest/lua_test.cc`
\n\n
If the folder where the file resides is not covered by
a `CMakeLists.txt` recursion then follow
the instructions for \ref DevManAddingCode "adding a code to the system".
\n\n
<B>Restrictions:</B> The lua-wrapper `.cc` file must be in a folder that
is not shared with `.cc` files that contain regular non-lua code. In the example
above it means we cannot add the wrapper `.cc` to the general `CHI_TECH` folder.
<I>This is because we use a script to convert the doxygen comments to
 text used in the input manual.</I>

## Step 2: Define a standard lua template with your input function name
In our case the test function will be `LuaTest`. We add `chi` to the front of
all input language to distinguish it from other lua functions.
\n\n
Lua-wrapper functions have no real arguments other than the lua-state and
always return an integer representing the number of items returned. For more
information on how to handle parameters in lua C-api consult the
<a href="https://www.lua.org/pil/24.html">lua-documentation</a>.
\n\n
<B>General Requirements:</B>
 - The doxygen comments must start with a forward dash
 and two asterisks. Like / * * but without spaces which was necessary to show
 it in this document.
 - If you use the `param` tag then it must be on its own line. The variable
 name is followed by an indication of the type of the variable and then
 by the description of the parameter. <I>This gets used in the documentation
 script to generate the input reference</I>.
 - The closing comment must be an asterisk forward slash. Like * / but without
 the space which was necessary to show it in this document. And must be on its
 own line.

\include "../../CHI_TECH/LuaTest/lua_test.cc"

## Step 3: Add it to the lua-register
In the project directory go to the file `CHI_TECH/ChiLua/chi_lua_register.h`
and add a `RegisterFunction` call to the register, registering the lua function
call to be connected with the lua console. <I>Without out this call the
wrapper will not be callable from lua.</I>
\n\n
If you would like for the function call to be listed under a category then
proceed the `RegisterFunction` call with a comment starting exactly with
`//module:` which will create a module with text matching the text following
this piece.
\n
\n
The module tabulation in the manual and the documentation for the parameters
will be processed with a script therefore it is important to comment
parameters and return values. Some specific notes in the lua-register file

 - `//module:` provides a categorization of the function on the documenation
 main page.
 - `//\ref` will provide a link anchor on the documenation main page.

\code
 .
 .
 .
RegisterFunction(chiLBSGroupsetSetGroupSubsets)
RegisterFunction(chiLBSGroupsetSetIterativeMethod)
RegisterFunction(chiLBSGroupsetSetResidualTolerance)
RegisterFunction(chiLBSGroupsetSetMaxIterations)
RegisterFunction(chiLBSGroupsetSetGMRESRestartIntvl)
RegisterFunction(chiLBSGroupsetSetWGDSA)
RegisterFunction(chiLBSGroupsetSetTGDSA)

//module:Test scripts
RegisterFunction(chiLuaTest)


\endcode

## Step 4: Include it in the documentation script
The module inclusion script is at `CHI_TECH/ChiLua/module_lua_inclusion`. Add
the folder where your lua-wrapper resides to this script:

\code
 .
 .
 .
Add_Folder("../../Modules/MonteCarlon/lua")
Add_Folder("../../Modules/DiffusionSolver/lua")
Add_Folder("../../Modules/LBSSteadyState/lua")

Add_Folder("../LuaTest")
\endcode


## Final step: Test it!!

Create a script input file and call your function. If it is linked incorrectly
the lua scripting system will throw an error like

`Test.lua:2: attempt to call a nil value (global 'chiLuaTest')`
 */
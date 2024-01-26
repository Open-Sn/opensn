# User Tutorials

## Lua for users

OpenSn input files are written in Lua: 
- Some of the code objects are directly accessible in the Lua console,
- Portion of the input can often be generated using Lua scripting.

Hence, we recommend the users to familiarize themselves with Lua. 
See the [official lua page](https://www.lua.org/) for more details.

1. A table of all existing lua functions accessible from the lua input file 
can be found on the main page of the [documentation](https://xxx.io/index.html), 
under the **Input Functions** label;
2. If you prefer a hierarchical list of those lua functions, there are also found from the main page of the 
documentation, by clicking on the left column on [
Input Reference](https://xxx.io/modules.html).

If you build the documentation locally, the above instructions should 
help you locate the local versions of those pages.

## Meshing

OpenSn handles many types of meshes, including arbitrary polyhedral meshes. 
Complex meshes are often generated externally and loaded into OpenSn. 
Hence, it is important for users to be familiar with CAD and mesh generators.

We provide some [meshing tutorial here](https://xxx.io/d4/d25/_mesh_tutorials.html) 
but recommend users seek additional knowledge when needed.

## Example of lua input files
A set of lua example files can be found:
1. in the documentation: In the left column on the main documentation page, click
on [Tutorials](https://xxx.io/d8/d2d/_tutorial00.html); 
these are fully commented out input files. 
2. in the ```tests``` folder: these files are the regression suite input files but 
they cover quite a few situations.

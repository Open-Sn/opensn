#!/bin/bash
# ================================ Check we are in the doc folder
cwd_name=${PWD##*/}          # to assign to a variable
cwd_name=${cwd_name:-/}        # to correct for the case where PWD=/

if [ $cwd_name != "doc" ]; then
  echo "ERROR: Script must be run from \"doc\" directory"
  exit 1
fi

# ================================ Check that we can operate the generated_files
#                                  folder
[ ! -d "generated_files/" ] && mkdir generated_files
cd generated_files/

cwd_name=${PWD##*/}          # to assign to a variable
cwd_name=${cwd_name:-/}        # to correct for the case where PWD=/

if [ $cwd_name != "generated_files" ]; then
  echo "ERROR: Director \"generated_files\" could not be made"
  exit 1
else
  rm -rf *
fi

# ================================ Generated static register documentation
cd ../../
if [ ! -f ./bin/ChiTech ]; then
  echo "ChiTech executable not found in bin folder"
  echo "Cannot update static registered items"
else
  echo "Running ChiTech"
  ./bin/ChiTech --suppress_color --dump-object-registry -v 2 > doc/generated_files/test.txt
  python3 "doc/scripts/MakeInputParamsDocs.py"
fi

#============================== Transform lua wrappers for documentation
# Lua wrapper functions are normally int chiFunction(lua_State* L)
# The script below uses their doc-strings to transform to doxy style
# functions
python3 "doc/scripts/BuildListOfLuaWrappers.py"

#============================== Build mainpage quick reference
python3 "doc/scripts/MakeMainPage.py"

#============================== Making main documentation
if [ "$1" == "input_doc_only" ]; then
  echo "Doing input documentation"
  grep -v "INPUT " doc/DoxyfileLua > doc/DoxyfileLuaLean
  echo 'INPUT = "doc/generated_files" "doc/PAGES"' >> doc/DoxyfileLuaLean
  for dee in $(find modules -type d -name "doc");
  do echo "INPUT += $dee" >> doc/DoxyfileLuaLean; done
  for dee in $(find framework -type d -name "doc");
  do echo "INPUT += $dee" >> doc/DoxyfileLuaLean; done
  echo "Running DoxyGen"
  doxygen "doc/DoxyfileLuaLean"
  rm doc/DoxyfileLuaLean
else
  echo "To generate only input documentation add the argument input_doc_only"
  echo "Example: ./YReGenerateDocumentation.sh input_doc_only"
  echo "Running DoxyGen"
  doxygen "doc/DoxyfileLua"
fi


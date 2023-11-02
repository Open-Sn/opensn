"""
This script walks the modules- and framework-folders and
looks for directories named "lua" or "LuaTest". When found,
it adds all the .cc files in those folders to a list.
"""
import os
import re

if not os.path.exists("doc/generated_files"):
    os.mkdir("doc/generated_files")

folders_to_walk = ["framework", "modules"]

wrapper_sources = []
for search_folder in folders_to_walk:
    for dir_path, sub_dirs, files in os.walk(search_folder):
        dir_name = os.path.basename(dir_path)
        if dir_name == "lua" or dir_name == "LuaTest":
            for file_name in files:
                base_name, extension = os.path.splitext(file_name)
                if extension == ".cc":
                    wrapper_sources.append(dir_path + "/" + file_name)

input_file_folders = ["test", "doc/examples"]
"""
The following function finds all lua input files in the 
input file folders. It then creates proxies for these files.
"""
input_files = []
main_input_proxies_file = open("doc/generated_files/input_wrappers_main.h", "w")
main_input_proxies_file.write("/** \page LuaInputExamples Lua Input Examples \n\n")

for input_file_folder in input_file_folders:
    # r=root, d=directories, f=files
    for r, d, f in os.walk(input_file_folder):
        for file_name in f:
            if '.lua' in file_name:
                file_path = r + "/" + file_name
                input_files.append(file_path)

                proxy_name = file_path.replace("/", "_")
                proxy_name = proxy_name.replace(".", "_")

                main_input_proxies_file.write("\subpage " + proxy_name + "  \n")

                dirs = r.split("/")
                incremental_root = "doc/generated_files"
                for dir in dirs:
                    incremental_root += "/" + dir
                    if not os.path.exists(incremental_root):
                        os.mkdir(incremental_root)

                input_proxies_file = open("doc/generated_files/" + r + "/" + \
                                          proxy_name + ".h", "w")

                input_proxies_file.write("/** \page " + proxy_name + " " + \
                                         file_path + "\n")
                input_proxies_file.write("\ingroup LuaInputExamples\n\n")
                input_proxies_file.write("\code\n")

                temp_file = open(file_path, "r")
                lines = temp_file.readlines()
                for line in lines:
                    input_proxies_file.write(line)
                temp_file.close()

                input_proxies_file.write("\endcode\n")
                input_proxies_file.write("*/\n\n\n")

                input_proxies_file.close()

main_input_proxies_file.write("*/\n\n")
main_input_proxies_file.close()

"""
The following function searches for the use of a function name in all
the lua input files
"""


def SearchForFunction(function_name, input_files):
    file_list = []

    for input_file in input_files:
        file = open(input_file, "r")
        content = file.read()

        if function_name in content:
            file_list.append(input_file)

        file.close()

    return file_list


"""
Next each file is scanned. Doxygen comments are extracted as well
as functions with the signature `int chiFunction(lua_State* L)`. The 
latter is identified with 3 items, i.e., the `int`, the prefix `chi` and
the presence of `lua_State`.

The functions only (without the doxy comments) are wrapped in a
namespace (chi_lua) and written to 
`doc/generated_files/lua_namespace.h` (the declaration file).

The doxy comments AND the functions are written to 
`doc/generated_files/lua_functions.c` (the definition file). 

For the definitions the doxy comments of a function is scanned for the presence 
of `param` and `return`. If `param` is encountered, as word0, then word1 is 
taken as the parameter name and word2 is taken as the parameter type. 
If `return` is encountered, as word0, then word1 is taken as the return type.

The function (declaration and definition) is then reformatted when written, such
that the return type replaces the classical `int` return and the parameters
replace the `lua_State* L` parameter.
"""
definition_file = open("doc/generated_files/lua_functions.c", "w")
declaration_file = open("doc/generated_files/lua_namespace.h", "w")
definition_file.write("namespace chi_lua\n{\n")
declaration_file.write("namespace chi_lua\n{\n")
for src_path in wrapper_sources:
    file = open(src_path, "r")

    definition_file.write("//" + src_path + "\n")
    lines = file.readlines()

    in_comment = False
    params = []
    returns = []
    comment_lines = []
    function_line = ""
    for line in lines:

        if line.find("/**") >= 0:
            in_comment = True

        if in_comment:
            comment_lines.append(line)

        if in_comment:
            param_start = line.find(r"\param")
            if param_start >= 0:
                words = line[param_start:].split()
                if len(words) >= 3:
                    param_type = words[2]
                    param_name = words[1]
                    if param_type[-1] == '.':
                        param_type = param_type[:(len(params)-2)]
                    params.append(param_type + " " + param_name)
            return_start = line.find(r"\return")
            if return_start >= 0:
                words = line[return_start:].split()
                if len(words) >= 2:
                    returns.append(words[1])

        if line.find("*/") >= 0:
            in_comment = False

        int_find = line.find("int")
        chi_find = line.find("chi")
        state_find = line.find("lua_State")

        if int_find >= 0 and chi_find >= 0 and state_find >= 0 and \
                (state_find > chi_find > int_find):
            words = re.split(r"\(|\)|\s", line.strip())

            # Build function name and filelist
            function_name = words[1]
            file_list = SearchForFunction(function_name, input_files)


            def WriteUsageExamples():
                if len(file_list) > 0:
                    definition_file.write("### Usage Examples:\n")
                    for usage_file in file_list:
                        proxy_name = usage_file.replace("/", "_")
                        proxy_name = proxy_name.replace(".", "_")
                        definition_file.write("\\ref " + proxy_name + "  \n")


            if len(returns) > 0:
                function_line = returns[0] + " "
            else:
                function_line = "void "

            function_line += function_name + "("
            for k in range(0, len(params)):
                function_line += params[k]
                if k < (len(params) - 1):
                    function_line += ", "
            function_line += ")"

            num_comment_lines = len(comment_lines)
            if num_comment_lines == 1:
                line = comment_lines[0].replace("*/", "")
                definition_file.write(line)
                WriteUsageExamples()
                definition_file.write("*/\n")
            elif num_comment_lines > 1:
                for k in range(0, num_comment_lines - 1):
                    definition_file.write(comment_lines[k])
                WriteUsageExamples()
                definition_file.write(comment_lines[num_comment_lines - 1])

            definition_file.write(function_line + "{}\n")
            declaration_file.write(function_line + ";\n")

            params = []
            returns = []
            comment_lines = []
            function_line = ""

    file.close()

definition_file.write("}//namespace chi_lua\n")
definition_file.close()

declaration_file.write("}//namespace chi_lua\n")
declaration_file.close()

import os
import textwrap

print("Creating InputParamsDocs")
main = open("doc/generated_files/test.txt", "r")

obj_list = []
object_dict = {"optional_params": [], "required_params": []}
param_dict = {}


def EverythingAfter(line_string: str, word: str):
    start = line_string.find(word)
    end = start + len(word)

    return line_string[end:]


lines = main.readlines()
for ell in range(0, len(lines)):
    line = lines[ell]
    words = line.split()

    if len(words) <= 1:
        continue

    if words[1] == "OBJECT_BEGIN":
        object_dict["name"] = words[2]
        object_dict["constructible"] = True

    if words[1] == "CLASS_NAME":
        if len(words) >= 3:
            end_of_word = line.find("CLASS_NAME") + len("CLASS_NAME")
            object_dict["class_name"] = line[end_of_word:]

    if words[1] == "DOC_GROUP":
        num_doc_groups = len(words) - 2
        object_dict["doc_groups"] = []
        for w in range(0, num_doc_groups):
            object_dict["doc_groups"].append(words[2 + w])

    if words[1] == "LUA_FUNCWRAPPER_BEGIN":
        object_dict["name"] = words[2]
        object_dict["function_wrapper"] = True

    if words[1] == "NOT_CONSTRUCTIBLE":
        object_dict["constructible"] = False

    if words[1] == "SYNTAX_BLOCK":
        object_dict["function_wrapper"] = False

    if words[1] == "OBJECT_END":
        obj_list.append(object_dict)
        object_dict = {"optional_params": [], "required_params": []}

    if words[1] == "LUA_FUNCWRAPPER_END":
        obj_list.append(object_dict)
        object_dict = {"optional_params": [], "required_params": []}

    if words[1] == "DESCRIPTION_BEGIN":
        object_dict["description"] = ""
        end_reached = False
        k = ell + 1
        while not end_reached:
            if lines[k].find("DESCRIPTION_END") >= 0:
                end_reached = True
                break  # from while

            object_dict["description"] += lines[k]
            k += 1

    if words[1] == "PARAM_BEGIN":
        param_dict["name"] = words[2]

    if words[1] == "PARAM_END":
        if param_dict["TAG"] == "OPTIONAL":
            object_dict["optional_params"].append(param_dict)
        else:
            object_dict["required_params"].append(param_dict)
        param_dict = {}

    if words[1] == "TYPE":
        param_dict[words[1]] = words[2]

    if words[1] == "TAG":
        param_dict[words[1]] = words[2]

    if words[1] == "DOC_STRING_BEGIN":
        param_dict["DOC_STRING"] = ""
        end_reached = False
        k = ell + 1
        while not end_reached:
            if lines[k].find("DOC_STRING_END") >= 0:
                end_reached = True
                break  # from while

            param_dict["DOC_STRING"] += lines[k]
            k += 1

    if words[1] == "DEFAULT_VALUE":
        param_dict[words[1]] = EverythingAfter(line, words[1]).strip()

    if words[1] == "CONSTRAINTS":
        param_dict[words[1]] = EverythingAfter(line, words[1]).strip()

    if words[1] == "LINKS":
        if len(words) > 2:
            param_dict[words[1]] = words[2:]

    ell += 1

main.close()

if not os.path.exists("doc/generated_files"):
    os.mkdir("doc/generated_files")

if not os.path.exists("doc/generated_files/input_docs"):
    os.mkdir("doc/generated_files/input_docs")

input_file_folders = ["test", "doc/examples"]
"""
The following function finds all lua input files in the 
input file folders.
"""
input_files = []
for input_file_folder in input_file_folders:
    # r=root, d=directories, f=files
    for r, d, f in os.walk(input_file_folder):
        for file_name in f:
            if '.lua' in file_name:
                file_path = r + "/" + file_name
                input_files.append(file_path)


def SearchForFunction(function_name, input_files):
    """
    The following function searches for the use of a function name in all
    the lua input files
    """
    file_list = []

    for input_file in input_files:
        file = open(input_file, "r")
        content = file.read()

        if function_name in content:
            file_list.append(input_file)

        file.close()

    return file_list


# The script below is a javascript that allows the collapsible boxes
# of parameters to function properly
dropdn_list_script2 = '''
{
    var coll = document.getElementsByClassName("droppy2");
    var i;
    
    for (i = 0; i < coll.length; i++) {
      coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = this.nextElementSibling;
        if (content.style.display === "block") {
          content.style.display = "none";
          this.childNodes[0].innerText = "►"
        } else {
          content.style.display = "block";
          this.children[0].innerText = "▼"
        }
      });
    }
}'''

# background-color: #313335;
# width: 100%;
# padding: 2px 2px;

button_style = '''
"
background-color: #edf0f5;
color: #9373A5;
cursor: pointer;
border: none;
text-align: left;
outline: none;
font-size: 15px;
margin: 2px 2px;
border: 2px solid #687372;
display: block;
width: 100%;
vertical-align: middle;
padding-top: 4px;
padding-bottom: 0px;
"
'''

sub_doc_button_style = '''
"
background-color: #ffffff;
color: #000000;
cursor: pointer;
border: none;
text-align: left;
outline: none;
font-size: 15px;
margin: 2px 2px;

display: block;
vertical-align: middle;
padding-top: 4px;
padding-bottom: 0px;
"
'''


# <span class="arrow" style="padding-left: 0px;">▼
# <span class="arrow" style="padding-left: 16px;">►</span>


def WriteHTMLParameters(file, obj_name, params, padding=0):
    """ This function writes HTML formatted parameters with
    dropdown boxes."""

    file.write('<div style="display: block;">\n')
    droppy = "droppy2"
    for param in params:
        doc_string = param["DOC_STRING"]

        # Build a brief doc string (Just the first sentence)
        doc_brief_end = doc_string.find(".")
        if doc_brief_end < 0:
            doc_str_end = len(doc_string)
        doc_brief = doc_string[0:doc_brief_end]

        # Filter out documentation for sub objects $(<object_name>)$
        sub_doc_begin = doc_string.find("$(")
        sub_doc_end = doc_string.find("$)")

        sub_doc_name = ""
        # if sub_doc_begin >= 0 and sub_doc_end >= 0:
        #     sub_string = doc_string[sub_doc_begin:(sub_doc_end + 2)]
        #     doc_string = doc_string.replace(sub_string, "")
        #
        #     sub_doc_name = sub_string[2:-2]
        if sub_doc_begin >= 0 and sub_doc_end >= 0:
            print(f'WARNING[{obj_name}]: Deprecated use of "$(" or ")$" in '
                  f'parameter \"{param["name"]}\"')

        if "LINKS" in param:
            if len(param["LINKS"]) > 0:
                sub_doc_name = param["LINKS"][0]

        file.write('<button type="button" class="' + droppy + '" ' +
                   'style=' + button_style + '>' +
                   '<span class="arrow" style="padding-left: 0px;">►' +
                   '</span><TT><B>' +
                   param["name"] + '</B></TT><span style="color: #c4c1c0;">' +
                   "&nbsp &nbsp &nbsp &nbsp" + doc_brief + '</span></button>\n')

        file.write('<div class="content" style="display: none;">\n')
        file.write('  <p><I>type=</I><span style="color: blue;"><TT>' +
                   param["TYPE"] + '</TT></span>. ')
        file.write(doc_string + "</p>\n")
        if "DEFAULT_VALUE" in param:
            file.write('  <p>Default value: <TT style="color:grey">' +
                       param["DEFAULT_VALUE"] +
                       "</TT></p>\n")
        if "CONSTRAINTS" in param:
            file.write("  <p>Allowable values: <TT>" + param["CONSTRAINTS"] +
                       "</TT></p>\n")

        if len(sub_doc_name) > 0:
            sub_obj = {}
            for obj in obj_list:
                if obj["name"] == sub_doc_name:
                    sub_obj = obj

            if len(sub_obj) == 0:
                continue

            pad = padding + 20
            file.write(f'<div class="content" style="padding-left: {pad}px;">\n')
            file.write('<HR>')
            if len(sub_obj["required_params"]) > 0:
                file.write('<button type="button" class="' + droppy + '" ' +
                           'style=' + sub_doc_button_style + '>' +
                           '<span class="arrow" style="padding-left: 0px;">►' +
                           '</span><B>' +
                           f"Required parameters for {sub_doc_name}" +
                           '</B></button>\n')
                file.write('<div class="content" style="display: none;'
                           'padding-left: 20px;">\n')
                WriteHTMLParameters(file, sub_doc_name, sub_obj["required_params"], pad)
                file.write('</div>\n')
            if len(sub_obj["optional_params"]) > 0:
                file.write('<button type="button" class="' + droppy + '" ' +
                           'style=' + sub_doc_button_style + '>' +
                           '<span class="arrow" style="padding-left: 0px;">►' +
                           '</span><B>' +
                           f"Optional parameters for {sub_doc_name}" +
                           '</B></button>\n')
                file.write('<div class="content" style="display: none;'
                           'padding-left: 20px;">\n')
                WriteHTMLParameters(file, sub_doc_name, sub_obj["optional_params"], pad)
                file.write('</div>\n')
            file.write('</div>\n')

        file.write('</div>\n\n')

    file.write('</div>\n')


for obj in obj_list:
    num_args = len(obj["required_params"]) + len(obj["optional_params"])
    reduced_name = obj["name"].replace(":", "_")
    file = open("doc/generated_files/input_docs/" + reduced_name + ".h", "w")
    dotted_name = reduced_name.replace("__", ".")

    class_name = dotted_name
    if "class_name" in obj:
        class_name = obj["class_name"]

    file.write("/** \\defgroup " + reduced_name + " " + class_name + "\n\n\n")

    if "doc_groups" in obj:
        doc_groups = obj["doc_groups"]
        if len(doc_groups) == 0:
            file.write("\\ingroup DocExperimental\n")
        else:
            for g in range(0, len(doc_groups)):
                file.write("\\ingroup " + doc_groups[g] + "\n")

    if "description" in obj:
        file.write(obj["description"] + "\n")
        if obj["description"].find("\ingroup") >= 0:
            print(f'WARNING[{reduced_name}]: instead of using \"\\ingroup in '
                  f'description use the SetDocGroup method')
        if obj["description"].find("\defgroup") >= 0:
            print(f'WARNING[{reduced_name}]: do not use \"\\defgroup in '
                  f'description. Rather use the method SetClassName')

    if "constructible" in obj:
        if obj["constructible"]:
            file.write(textwrap.dedent('''
            ## Example usage:
            Create this object:
            \\code
            params = 
            {
                param_name1 = value1,
                param_name2 = value2,
                --etc.
            }\n'''))
            file.write(obj["name"].replace("::", ".") + ".Create(params)\n")
            file.write("\\endcode\n")
        else:
            file.write("<B>Note:</B> This object is not constructable\n")

    if "function_wrapper" in obj:
        if obj["function_wrapper"]:
            file.write(textwrap.dedent('''
                ## Function Syntax:
                \\code
                '''))
            file.write(obj["name"].replace("::", ".") + "(")
            for p in range(0, num_args):
                type_str = "UNKNOWN_TYPE"
                for param in obj["required_params"]:
                    if param["name"] == f"arg{p}":
                        type_str = param["TYPE"]
                        break
                for param in obj["optional_params"]:
                    if param["name"] == f"arg{p}":
                        type_str = param["TYPE"]
                        break

                file.write(f"{type_str} arg{p}")
                if p < num_args - 1:
                    file.write(", ")
            file.write(")\n")
            file.write("\\endcode\n")

    if len(obj["required_params"]) > 0:
        file.write("## Required Input parameters\n")
        file.write('\\htmlonly\n\n')
        WriteHTMLParameters(file, dotted_name, obj["required_params"])
        file.write('\\endhtmlonly\n\n')

    if len(obj["optional_params"]) > 0:
        file.write("## Optional Input parameters\n")
        file.write('\\htmlonly\n\n')
        WriteHTMLParameters(file, dotted_name, obj["optional_params"])
        file.write('\\endhtmlonly\n\n')

    file.write('\\htmlonly\n\n')
    file.write('<script>\n')
    file.write(dropdn_list_script2)
    file.write('</script>\n\n')
    file.write('\\endhtmlonly\n\n')

    file_list = SearchForFunction(obj["name"].replace("::", "."), input_files)
    if len(file_list) > 0:
        file.write("### Usage Examples:\n")
        for usage_file in file_list:
            proxy_name = usage_file.replace("/", "_")
            proxy_name = proxy_name.replace(".", "_")
            file.write("\\ref " + proxy_name + "  \n")

    file.write("*/")
    file.close()

print("Done creating InputParamsDocs")

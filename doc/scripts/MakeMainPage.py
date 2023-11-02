"""
This script reads the mp_header.txt and mp_tableinput.txt files.
It then first echos mp_header.txt directly to `doc/generated_files/MainPage.h`.
Thereafter it processes mp_tableinput.txt line by line.
- If the word `module:` is encountered then the text following that is printed
  as a level 3 heading (i.e. ###)
- If the word `print:` is encountered then the text following that is printed
  as is.
- If the word `function:` is encountered then the function is collected into
  a list that is formatted as a table.
- If the word `include:` is encountered then the file pointed to by the path
  following this word, will be processed in the same fashion
"""

output_file = open("doc/generated_files/MainPage.h", "w")

# ==================================================== Read & echo header
file = open("doc/scripts/mainpage_header.txt", "r")

lines = file.readlines()
for line in lines:
    output_file.write(line)

file.close()

# ==================================================== Process quickref tables
lines = []


def RecursiveIncludes(file_path):
    file_obj = open(file_path, "r")
    file_lines = file_obj.readlines()
    file_obj.close()

    for line in file_lines:
        lines.append(line)
        words = line.split()
        if len(words) > 0:
            first_word = words[0]
            if first_word == "include:":
                file_path = line[len(first_word):].strip()
                RecursiveIncludes(file_path)


RecursiveIncludes("doc/scripts/mainpage_quickref.txt")

# ===================================== Function for writing a table
def write_table(table_entries):
    N = len(table_entries)
    if N == 0:
        return
    output_file.write("<table>\n")

    N_mod_3 = N % 3

    for n in range(0, N-N_mod_3, 3):
        output_file.write("<tr>\n")
        for k in range(0, 3):
            output_file.write('<td width="33%">')
            output_file.write(table_entries[n + k])
            output_file.write("</td>\n")
        output_file.write("</tr>\n")

    output_file.write("<tr>\n")
    for n in range(N-N_mod_3, N):
        output_file.write('<td width="33%">')
        output_file.write(table_entries[n])
        output_file.write("</td>\n")

    output_file.write("</tr>\n")

    output_file.write("</table>\n\n")

table_entries = []
section_num = 1
sub_section_num = 1
sub_sub_section_num = 1
for line in lines:

    words = line.split()
    if len(words) == 0:
        continue

    first_word = words[0]

    if line.find("\section") >= 0:
        section_num += 1
        sub_section_num = 1

    if first_word == "module:":
        output_file.write(
            f"\\subsection MainPage{section_num}_{sub_section_num} " +
            f"{section_num}.{sub_section_num} " +
            line[len(first_word):].strip() + "\n")
        sub_section_num += 1

    if first_word == "submodule:":
        output_file.write(
            f"\\subsubsection MainPage{section_num}_" +
            f"{sub_section_num}_{sub_sub_section_num} " +
            line[len(first_word):].strip() + "\n")
        sub_sub_section_num += 1

    if first_word == "print:":
        output_file.write(line[len(first_word):].strip() + "\n\n")

    if first_word == "function:":
        table_entries.append(line[len(first_word):].strip()+"()")

    if first_word == "module_end":
        write_table(table_entries)
        table_entries = []


output_file.write("\n*/")
output_file.close()

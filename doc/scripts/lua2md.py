# Convert .lua file with in-code documentation into a .md file

import argparse
import os.path
import re

# Content of the generated markdown file
markdown = ""
# Lua source that will be included at the end of the file
lua_src = ""


t_DOCO = r'(?sm)--\[\[(.+?)--\s*\]\]'
t_CODE = r'.*?\n'


def extract_comment_content(str):
    """
    Extract content from a the in-code doco comment `--[[ @doc text --]]`
    """
    m = re.search(r'(?sm)--\[\[\s+@doc(.+)--\s*\]\]', str)
    if m:
        return m.group(1).strip()
    else:
        return ""


def remove_token(s, token):
    n = len(token)
    return s[n:]


def parse(src):
    global markdown
    global lua_src

    # state variable for indicating if we are inside a code block or not
    code_block = False

    s = src
    while len(s) > 0:
        if re.match(t_DOCO, s):
            # close opened code block
            if code_block:
                markdown += "```\n"
                markdown += "\n"
                code_block = False

            m = re.match(t_DOCO, s)
            token = m.group(0)
            s = remove_token(s, token)

            markdown += extract_comment_content(token) + "\n"
            markdown += "\n"
        elif re.match(t_CODE, s):
            m = re.match(t_CODE, s)
            token = m.group(0)
            s = remove_token(s, token)

            code = token.rstrip()
            if len(code) > 0:
                # open a code block
                if not code_block:
                    markdown += "```lua\n"
                    code_block = True

                markdown += code + "\n"

            lua_src += token
        else:
            raise SystemExit("Nothing matched")

    # if last block was a code block, close it
    if code_block:
        markdown += "```\n"
        markdown += "\n"


# main
parser = argparse.ArgumentParser(prog='lua2md')
parser.add_argument('lua_filename')
parser.add_argument('md_filename')
parser.add_argument('-d', '--root-directory', dest='root_directory', default=None,
                    help='root directory')
args = parser.parse_args()

with open(args.lua_filename, "r") as f:
    src = ''.join(f.readlines())
    parse(src)

if args.root_directory is not None:
    display_file_name = os.path.relpath(args.lua_filename, args.root_directory)
else:
    display_file_name = args.lua_filename
with open(args.md_filename, "w") as f:
    f.write(markdown)
    # append the complete input file
    f.write("___\n")
    f.write("## The complete input is below:\n")
    f.write(
        "You can copy/paste the text below or look in the file named ```{}```:\n".format(
            display_file_name))
    f.write("```lua\n")
    f.write(lua_src)
    f.write("```\n")

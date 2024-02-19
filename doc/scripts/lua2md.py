# Convert .lua file with in-code documentation into a .md file

import argparse
import re
import ply.lex as lex
import ply.yacc as yacc

# Content of the generated markdown file
markdown = ""
# Lua source that will be included at the end of the file
lua_src = ""


def extract_comment_content(str):
    """
    Extract content from a the in-code doco comment `--[[ @doc text --]]`
    """
    m = re.search(r'(?sm)--\[\[\s+@doc(.+)--\s*\]\]', str)
    if m:
        return m.group(1).strip()
    else:
        return ""


# lexer
tokens = [
    'DOCO',
    'CODE'
]

t_DOCO = r'(?sm)--\[\[(.+?)--\s*\]\]'
t_CODE = r'.*?\n'


def t_error(t):
    t.lexer.skip(1)


def p_file(p):
    'file : blk'
    pass


def p_file_seq(p):
    'file : file blk'
    pass


def p_doco(p):
    'blk : DOCO'
    global markdown

    markdown += extract_comment_content(p[1]) + "\n"
    markdown += "\n"


def p_code(p):
    'blk : code_block'
    global markdown
    global lua_src

    code = p[1].strip()
    if len(code) > 0:
        markdown += "```\n"
        markdown += code + "\n"
        markdown += "```\n"
        markdown += "\n"

    lua_src += p[1]


def p_code_block_one(p):
    'code_block : CODE'
    p[0] = p[1]
    if p[1] != '\n':
        p[0] = p[1]
    else:
        p[0] = ""


def p_code_block_seq(p):
    'code_block : code_block CODE'
    p[0] = p[1] + p[2]


# main
parser = argparse.ArgumentParser(prog='lua2md')
parser.add_argument('lua_filename')
parser.add_argument('md_filename')
args = parser.parse_args()

with open(args.lua_filename, "r") as f:
    src = ''.join(f.readlines())
    lexer = lex.lex()
    parser = yacc.yacc()
    result = parser.parse(src)

with open(args.md_filename, "w") as f:
    f.write(markdown)
    # append the complete input file
    f.write("___\n")
    f.write("## The complete input is below:\n")
    f.write(
        "You can copy/paste the text below or look in the file named ```{}```:\n".format(
            args.lua_filename))
    f.write("```\n")
    f.write(lua_src)
    f.write("```\n")

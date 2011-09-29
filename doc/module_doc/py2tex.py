# py2tex.py -- Translate Python source code to % LaTeX
# \LaTeX\ code that can be typeset using the \|py2tex| documentstyle
# option.
# \ESC
# To typeset a Python module called \|foo.py| with py2tex, create a
# \LaTeX\ file along the following lines.
#
# \begin{verbatim}
# %% frame.tex -- wrapper around foo
# \documentstyle[...,py2tex,...]{...}
# ...
# \begin{document}
# ...
# \PythonSource{foo.pt}
# ...
# \end{document}
# \end{verbatim}
# Then give the command
#
# \I{}{1}{\tt\$ py2tex -o foo.pt foo.py}
#
# Finally run \LaTeX\ on the previously constructed wrapper, like this
#
# \I{}{1}{\tt\$ latex frame}
#
# This will give you a \|.dvi| file that you can print in the normal
# way.
#
# Note that normally the comments are interpreted by \LaTeX.  This
# allows for formulae and other fancy stuff.  However, if you don't
# need this, or if you want to typeset programs that were not
# specifically written to be typeset with py2tex, you can leave
# comments uninterpreted by calling the \|py2tex| script with
# the~\|-v| option.  The same effect can be obtained by ending a
# comment with `\|%ASCII|'. It is also possible to switch back to
# interpreted mode by inserting a comment ending in `\|%TeX|' or
# `\|%LaTeX|'.
#
# Here are some guidelines for writing Python code to be typeset using
# py2tex. Each line of Python code is typeset by \LaTeX\ as a paragraph
# where, in case it is broken up into more than one line, all lines
# following the first are indented by one and a half standard
# indentation more than the indentation of the first line. Py2tex does
# not count parentheses to determine whether a line is a continuation
# of the previous or not.  So if you want it to be indented
# appropriately, escape the end of the previous line with a backslash.
# Then py2tex will treat the joined lines as one line, and it will
# inform \LaTeX\ that the escaped line breaks are good points to break
# it up again.  Because \LaTeX\ may decide to break the code at other
# positions (or not at all), these lines will not be numbered.
#
# Consecutive lines that start with a single hash mark (\##) right
# after the indentation are joined and typeset in a \|\vbox| (more
# precise: a \|\vtop|).  This is called a block comment.  Indentation
# changes have no effect within a block comment.  It is possible to
# escape from the \|\vbox| and set the remainder of the block comment
# in what Knuth calls `outer vertical mode' by using the \|\ESC|
# command.  This can be used to incorporate long stretches of \LaTeX\
# code that can spread out over several pages.  Unindented block
# comments are automatically escaped in their entirety.
#
# If a line starts with at least two hash marks it is typeset as if it
# followed some Python code.  The second hash mark also switches
# immediately back to Python mode (see below).  This feature is also
# implemented for {\sc ascii} mode, while the general escape to Python
# mode is not.  (This feature is intended to disable lines of Python
# code by placing two hash marks before them.  This ensures that the
# formatting will be very similar to the uncommented version.)
#
# Comments following Python code are typeset on the same line as the
# Python code, separated from it by a \|\quad| space and the hash
# mark.
#
# Both in block and in line comments the hash mark is used to switch
# between \LaTeX\ and Python mode, just like the dollar sign (\$) is
# used to switch between horizontal and math mode.  This means
# that hash marks are not visible as such in the output.  However, two
# consecutive hash marks are passed to \LaTeX\ as one.  This means that
# it is possible to typeset a hash mark by putting \|\####| in a
# comment.  (This can also be used to define \LaTeX\ macros and to
# include \|\halign| templates, albeit at the expense of doubling all
# hash marks.)  Note that this works only in \LaTeX\ mode, {\em not\/}
# in {\sc ascii} mode.
#
# So if you type
#
# \vskip\parskip\vbox{\parskip=0mm
# \I{}{1}{\tt \## \% LaTeX}
# \I{}{1}{\tt \## Hash mark in comment:\ \|\####|,}
# \I{}{1}{\tt \## formula in comment:\ \$i\_0\char`\\to\char`\\infty\$.}
# \I{}{1}{\tt print chr (i) \## where \##040<=i<=0x7E}
# \I{}{1}{\tt \##\## print #'#'# \## print one hash sign \% ASCII}
# \I{}{1}{\tt \##\## print #i_0*'#'# \## where i\_0 is \##hash signs}
# \par}
#
# you get
#
# \vskip\parskip\vbox{\parskip=0mm
# \B{}{1}{
#    Hash mark in comment: \##,
#    formula in comment: $i_0\to\infty$.
# }
# \I{}{1}#print chr (i)#\quad\## where #040<=i<=0x7E#
# \I{}{1}\## #print '#'#\quad\## print one hash sign
# \I{}{1}\## #print i_0*'#'#\quad\## where i\_0 is \##hash signs
# \par}
#
# Triple quoted strings that occur as the first non-comment after a line
# that ends in a colon (:) are treated as documentation strings.
# There is are three different options for treating
# them. If #docprocess=='none'#, this results in the ``Same 'ol behaviour'':
# \PythonSource*{docstring-none.pt}
# If #docprocess=='plain'#, docstrings are typeset as verbatim
# comments except with thick solid lines instead of thin double
# lines:
# \PythonSource*{docstring-plain.pt}
# If #docprocess=='struct'#, docstrings are typeset as structured text
# as defined by the doc-sig. This is so people can potentially
# write programs that look good both under gendoc and py2tex.
# \PythonSource*{docstring-struct.pt}
#
# It is possible to include the formatted version of another Python
# source file using the \|\PythonSource*| macro.  This was done below
# to give an example of the use of class #Interpret#.
# The starred version of the macro is needed to drop the line numbers,
# otherwise they would be typeset through the lines that mark the
# block comment.  The starred version of \|\PythonSource| also drops
# the section heading. If you escape the block comment (using \|\ESC|)
# you can use the unstarred version again.
#
# \begin{table}\sf
# \hrule height 1mm
# \smallskip
# {\bf Table 1.} Some Python constructs get special typographic treatment
#
# \medskip
# \halign{\tt ##\hfil\qquad&##\hfil\cr
# \omit\bf Python\hfil        &\omit\bf \LaTeX \cr
# =                &#=# \cr
# ==                &#==# \cr
# <=, >=        &#<=, >=# \cr
# !=, <>        &#!=# \cr
# <<, >>        &#<<, >># \cr
# and, or, not        &#and, or, not# \cr
# in, not in        &#in, not in# \cr
# is, is not        &#is, is not# \cr
# }
# \smallskip
# \hrule height 1mm
# \end{table}
#
# Finally some remarks about the formatting of Python
# constructs. Identifiers (keywords, variables and functions) are
# typeset in sans serif.  If an identifier consists of only one
# character, it is typeset in {\em math italic\/} instead of sans
# serif.  Keywords are typeset in boldface, functions (actually:
# identifiers before opening parentheses) are typeset slanted.  These
# typefaces can be changed by redefining some of the macros in
# \|py2tex.sty|.  See the documentation of the style file for
# customization instructions.
#
# Some constructs that get special treatment are listed in Table~1.
# This special treatment is optional.  If the class is initialized
# with an extra argument that evaluates to false, or if the
# #no_math()# method is used, then no special treatment is done for
# these constructs.  (Special treatment can be turned back on half way
# through a file using the #math()# method.)
#
# In strings, characters outside the range #' '#--#'~'# are typeset as
# standard escape sequences (\eg, \|TAB| is typeset as #'        '#,
# \|ESC| is typeset as #''#).  A floating point literal with an
# exponent has its exponent written out as a power of ten (\eg,
# \|3e-6| is typeset as #3e-6#).  Hexadecimal literals are typeset in
# a typewriter font with a lower case \|x| and uppercase digits (\eg,
# \|0X007e| is typeset as #0X007e#).  Octal literals are typeset in
# italics (\eg, \|0377| is typeset as #0377#).

import os, re, string, sys, time

# Usage of class #Interpret#.
#\ESC
# \PythonSource*{example.pt}
# Note that #sys.stdin# is used if #name in (None, '-')#.
#
# The other methods can best be viewed as private to the class.
class Interpret:
    def __init__ (self, name, math = 1, interpret = 1, docprocess='none'):
        if name == None:
            self._name = '-'
        else:
            self._name = name
        if self._name == '-':
            self._name = '(stdin)'
            mtime = time.asctime (\
                  time.localtime (time.time ()))
            self._file = sys.stdin
        else:
            mtime = time.asctime (\
                  time.localtime (os.stat (name) [8]))
            self._file = open (self._name, 'r')
            self._name = os.path.basename (self._name)
        preamble = '\\File{%s}{%s}\n\n' % (self._name, mtime)
        if not math: preamble = preamble + '\\PythonNoMath\n\n'
        self._translation = [preamble,]
        self._docstuff = [preamble,]
        self._math = math
        self._line_nr = 0
        self._line = None
        self._old_line = None
        self._eof = 0
        self._indent_stack = [0]
        self._no_break = 0
        self._interpret_comments = interpret
        self._docprocess = docprocess
        self._docstring = 1
    def math (self):
        if not self._math:
            self._translation.append ('\\PythonMath\n')
            self._docstuff.append ('\\PythonMath\n')
            self._math = 1
    def no_math (self):
        if not self._math:
            self._docstuff.append ('\\PythonNoMath\n')
            self._math = 0
    def interpret (self):
        self._interpret_comments = 1
    def verbatim (self):
        self._interpret_comments = 0
    def close (self):
        self._file.close ()
        self._line_nr = 0
        self._line = None
        self._old_line = None
        self._indent_stack = [0]
        self._translation = []
        self._docstuff = []
        self._eof = 1
        self._no_break = 0
    def flush (self):
        self._file.flush ()
    def next_line (self):
        if self._old_line != None:
            self._line = self._old_line
            self._old_line = None
            self._line_nr = self._line_nr + 1
            return
        self._line = self._file.readline ()
        if self._line == '':
            self._eof = 1
            raise EOFError
        if self._line [-1] == '\n': self._line = self._line [:-1]
        self._line_nr = self._line_nr + 1
    def undo_line (self):
        if self._line != None:
            self._old_line = self._line
            self._line = None
            self._line_nr = self._line_nr - 1
    def close_tex (self, tex):
        while 1:
            if tex [-2:] == '\\ ':
                tex = tex [:-2]
            elif tex [-4:] == '\\BP ':
                tex = tex [:-4]
            else:
                break
        if tex not in ('$', '${}'):
            self._translation.append (tex + '$')
            if tex.find(r'\K{def}') != -1 or tex.find(r'\K{class}') != -1:
                self._docstuff.append(tex + '$')
            elif len(self._docstuff) > 0 and self._docstuff[-1].startswith('\\I{'):
                self._docstuff = self._docstuff[:-1]
    def tr_indentation (self):
        length = white_re.match (self._line)
        if length < 0: raise error
        indent = 0
        for c in self._line [:length]:
            indent = indent + 1
            if c == '\t':
                indent = indent + 8
                indent = indent & ~0x7
        self._line = self._line [length:]
        while indent < self._indent_stack [-1]:
            del self._indent_stack [-1]
        if indent > self._indent_stack [-1]:
            self._indent_stack.append (indent)
        self._indentation = len (self._indent_stack) - 1
    def tr_comment_line (self):
        if self._interpret_comments:
            length = verbatim_re.search (self._line)
            if length >= 0:
                self.verbatim ()
                self._line = self._line [:length]
            while 1:
                hash = string.find (self._line, '#')
                if hash >= 0:
                    if len (self._line) > hash + 1 and \
                          self._line [hash + 1] == '#':
                        self._translation.append (\
                              self._line [:hash] + '#')
                        #self._docstuff.append (\
                        #      self._line [:hash] + '#')
                        self._line = self._line [hash + 2:]
                        continue
                    self._translation.append (self._line [:hash])
                    #self._docstuff.append (self._line [:hash])
                    self._line = self._line [hash + 1:]
                    self.tr_code (0) # No continued lines in comments.
                    if len (self._line) <= 0: break
                    if self._line [0] != '#': raise error
                    self._line = self._line [1:]
                else:
                    break
            self._translation.append (self._line + '\n')
            #self._docstuff.append (self._line + '\n')
        else:
            length = interpret_re.search (self._line)
            if length >= 0:
                self.interpret ()
                self._line = self._line [:length]
            while len (self._line) > 0:
                length = ordinary_re.match (self._line)
                if length > 0:
                    self._translation.append (self._line [:length])
                    #self._docstuff.append (self._line [:length])
                if len (self._line) > length:
                    char = self._line [length]
                    if char in '<>\\{|}~':
                        self._translation.append (\
                              '{\\tt\\char`\\%s}' % char)
                        #self._docstuff.append (\
                        #      '{\\tt\\char`\\%s}' % char)
                    else:
                        self._translation.append ('\\' + char)
                        #self._docstuff.append ('\\' + char)
                self._line = self._line [length+1:]
            self._translation.append ('\n')
            #self._docstuff.append ('\n')
    def tr_block_comment (self):
        if self._line [0] != '#': raise error
        outer = self._indentation == 0
        if outer:
            if self._line_nr > 1:
                self._translation.append (\
                      '\\PythonOuterBlock\n')
                #self._docstuff.append (\
                #      '\\PythonOuterBlock\n')
            else:
                self._translation.append (\
                      '\\PythonOuterBlock*\n')
                #self._docstuff.append (\
                #      '\\PythonOuterBlock*\n')
        else:
            self._translation.append ('\\B{%d}{%d}{%%\n' % \
                  (self._line_nr, self._indentation))
            #self._docstuff.append ('\\B{%d}{%d}{%%\n' % \
            #      (self._line_nr, self._indentation))
        try:
            white = white_re.match (self._line, 1)
            if white < 0: raise error
            self._line = self._line [white:]
            while 1:
                self.tr_comment_line ()
                self.next_line ()
                white = white_re.match (self._line)
                if white < 0: raise error
                if len (self._line) > white and self._line [white] == '#' and \
                      self._line [white:white + 2] != '##':
                    self._line = self._line [white + 1:]
                    white = white_re.match (self._line)
                    if white > 0: self._line = self._line [white:]
                    continue
                self.undo_line ()
                return
        finally:
            if outer:
                self._translation.append (\
                      '\\PythonOuterBlockEnd\n')
                #self._docstuff.append (\
                #      '\\PythonOuterBlockEnd\n')
            else:
                self._translation.append ('}\n')
                #self._docstuff.append ('}\n')
    def tr_comment (self):
        self._translation.append ('\\#\\ ')
        #self._docstuff.append ('\\#\\ ')
        while self._line [:2] == '##':
            self._line = self._line [2:]
            self.tr_code (0) # No continued lines in comments.
            if self._line [:1] == '#':
                self._translation.append ('\\quad\\#\\ ')
                #self._docstuff.append ('\\quad\\#\\ ')
        if len (self._line) < 1:
            self._translation.append ('\n')
            #self._docstuff.append ('\n')
            return
        if self._line [0] != '#': raise error
        white = white_re.match (self._line, 1)
        if white < 0: raise error
        self._line = self._line [white:]
        self.tr_comment_line ()
        return
    def tr_string (self, token):
        quote = token [0]
        tl = len (token)
        self._translation.append ('\S{' + token)
        #self._docstuff.append ('\S{' + token)
        while 1:
            pos = string.find (self._line, quote)
            if pos > 0:
                self._translation.append ('\\verb*%s%s' % \
                      (quote, ctrl_protect (self._line [:pos + 1])))
                #self._docstuff.append ('\\verb*%s%s' % \
                #      (quote, ctrl_protect (self._line [:pos + 1])))
                if escape_re.match (self._line [:pos]) == pos:
                    self._translation.append (quote)
                    #self._docstuff.append (quote)
                    self._line = self._line [pos + 1:]
                    continue
                self._line = self._line [pos:]
                pos = 0
            if pos >= 0:
                if self._line [:tl] == token:
                    self._translation.append (token + '}')
                    #self._docstuff.append (token + '}')
                    self._line = self._line [tl:]
                    return
                self._translation.append (quote)
                #self._docstuff.append (quote)
                self._line = self._line [1:]
            else:
                self._translation.append ('\\verb*%s%s%s' % \
                      (quote, ctrl_protect (self._line), quote))
                #self._docstuff.append ('\\verb*%s%s%s' % \
                #      (quote, ctrl_protect (self._line), quote))
                self._line = ''
                if tl == 1:
                    self._translation.append ('}')
                    #self._docstuff.append ('}')
                    return
                self.next_line ()
                self._translation.append ('\n\\I{%d}{0}' % \
                      self._line_nr)
                #self._docstuff.append ('\n\\I{%d}{0}' % \
                #      self._line_nr)
        return

    def tr_docstring_plain (self):
        length = quote_re.match(self._line)
        token = self._line[:length]
        self._line = self._line[length:]
        quote = token [0]
        tl = len (token)

        if self._indentation == 0:
            self._translation.append ('\\PythonDocBlock\n')
            self._docstuff.append ('\\PythonDocBlock\n')
        else:
            self._translation.append ('\DS{%s}{%s}{%%\n' % 
                                      (self._line_nr, self._indentation))
            self._docstuff.append ('\DS{%s}{%s}{%%\n' % 
                                      (self._line_nr, self._indentation))
        while 1:
            pos = string.find (self._line, quote)
            if pos > 0:
                self._translation.append ('\\verb%s%s' % \
                      (quote, ctrl_protect (self._line [:pos+1])))
                self._docstuff.append ('\\verb%s%s' % \
                      (quote, ctrl_protect (self._line [:pos+1])))
                if escape_re.match (self._line [:pos]) == pos:
                    self._translation.append (quote)
                    self._docstuff.append (quote)
                    self._line = self._line [pos + 1:]
                    continue
                self._line = self._line [pos:]
                pos = 0
            if pos >= 0:
                if self._line [:tl] == token:
                    self._line = self._line [tl:]
                    break
                self._translation.append (quote)
                self._docstuff.append (quote)
                self._line = self._line [1:]
            else:
                self._translation.append ('\\verb%s%s%s' % \
                      (quote, ctrl_protect (self._line), quote))
                self._docstuff.append ('\\verb%s%s%s' % \
                      (quote, ctrl_protect (self._line), quote))
                self._line = ''
                if tl == 1:
                    break
                self.next_line ()
                # XXX This assumes 8 spaces per tab.
                wchars = white_re.match(self._line)
                spaces = re.sub('\t', ' '*8, self._line[:wchars])
                indent = white_re.match(spaces)
                ## print `spaces`, indent, self._indentation
                self._line = spaces[self._indentation*4:] + self._line[wchars:]
                self._translation.append ('\\\\ \n')
                self._docstuff.append ('\\\\ \n')
        if self._indentation == 0:
            self._translation.append ('\n\\PythonDocBlockEnd\n')
            self._docstuff.append ('\n\\PythonDocBlockEnd\n')
        else:
            self._translation.append ('}\n')
            self._docstuff.append ('}\\vskip 5pt \n')

    def tr_docstring_struct (self):
        length = quote_re.match(self._line)
        token = self._line[:length]
        self._line = self._line[length:]
        quote = token [0]
        tl = len (token)
        if self._indentation == 0:
            self._translation.append ('\\PythonDocBlock\n')
            self._docstuff.append ('\\PythonDocBlock\n')
        else:
            self._translation.append ('\DS{%s}{%s}{%%\n' % 
                                      (self._line_nr, self._indentation))
            self._docstuff.append ('\DS{%s}{%s}{%%\n' % 
                                      (self._line_nr, self._indentation))
        docstring = []
        while 1:
            pos = string.find (self._line, quote)
            if pos > 0:
                docstring.append (self._line [:pos])
                if escape_re.match (self._line [:pos]) == pos:
                    docstring.append (quote)
                    self._line = self._line [pos+1:]
                    continue
                self._line = self._line [pos:]
                pos = 0
            if pos >= 0:
                if self._line [:tl] == token:
                    self._line = self._line [tl:]
                    break
                docstring.append (quote)
                self._line = self._line [1:]
            else:
                docstring.append (self._line)
                self._line = ''
                if tl == 1:
                    break
                self.next_line ()
                docstring.append ('\n')
        docstring = string.joinfields(docstring, '')
        import struct2latex
        structstring = str(struct2latex.LaTeX(docstring))
        if self._indentation == 0:
            self._translation.append ('%s\n\\PythonDocBlockEnd\n' % structstring)
            self._docstuff.append ('%s\n\\PythonDocBlockEnd\n' % structstring)
        else:
            self._translation.append ('%s}' % structstring)            
            self._docstuff.append ('%s}' % structstring)            

    def tr_code (self, allow_continue = 1):
        tex = '$'
        try:
            careful   = 0
            while 1:
                white = white_re.match (self._line)
                if white > 0:
                    self._line = self._line [white:]
                if len (self._line) <= 0: return
                if self._line == '\\':
                    if allow_continue:
                        tex = tex + '\\BP '
                        self.next_line ()
                        continue
                    else:
                        self._line = ''
                        return
                if self._line [0] == '#': return
                length = token_re.match (self._line)
                if length < 1:
                    length = numeral_re.match (self._line)
                    if length < 1:
                        tex = tex + self._line [0]
                        self._line = self._line [1:]
                        careful = 0
                    else:
                        token = self._line [:length]
                        self._line = self._line [length:]
                        if careful: tex = tex + '\\ '
                        tex = tex + tr_numeral (token)
                        careful = 1
                    continue
                token = self._line [:length]
                self._line = self._line [length:]
                token = self.double (token)
                if token == ':':
                    self._docstring = 1
                else:
                    self._docstring = 0
                if token in ('{', '}'):
                    tex = tex + '\\' + token
                    careful = 0
                    continue
                if token in reserved_operators:
                    tex = tex + '\\O{%s}' % token
                    careful = 0
                    continue
                if token [0] in string.letters + '_':
                    if careful: tex = tex + '\\ '
                    new_careful = 1
                    if token in reserved:
                        if tex [-2:] not in ('$', '\\ ') and not careful:
                            tex = tex + '\\ '
                        tex = tex + '\\K{%s}' % token
                        if token == 'if': tex = tex + '\\,'
                        if token not in single: tex = tex + '\\ '
                        new_careful = 0
                    else:
                        token = usc_protect (token)
                        length = function_re.match (self._line)
                        if length > 0:
                            self._line = self._line [length:]
                            tex = tex + '\\F{%s}\\,(' % token
                            new_careful = 0
                        else:
                            if len (token) == 1:
                                tex = tex + token
                            else:
                                tex = tex + '\\V{%s}' % token
                    careful = new_careful
                    continue
                if token [0] in '\'"':
                    self.close_tex (tex + '{}')
                    self.tr_string (token)
                    tex = '${}'
                    careful = 0
                    continue
                if '{' in token or '}' in token:
                    raise ValueError, "brace in token '%s'" % token
                tex = tex + '\\Y{%s}' % token
                careful = 0
        finally:
            self.close_tex (tex)
    def double (self, token):
        if token not in ('not', 'is'): return token
        white = white_re.match (self._line)
        if white > 0:
            self._line = self._line [white:]
        next_length = token_re.match (self._line)
        if next_length > 0:
            next = self._line [:next_length]
            if (token, next) in (('not', 'in'), \
                  ('is', 'not')):
                self._line = self._line [next_length:]
                return token + ' ' + next
        return token
    # Method #translate()# is the interface to the #Interpret# class.
    # It calls the $\F{tr\_}xxx()$ methods to process indentation, code,
    # comments and strings.
    def translate (self):
        self._docstuff = []
        self._translation = []
        if self._eof: return None
        try:
            empty = 0
            self.next_line ()
            while white_re.match (self._line) == len (self._line):
                empty = empty + 1
                self.next_line ()
            if empty > 0:
                self._translation.append ('\\E{%d}' % empty)
                #self._docstuff.append ('\\E{%d}' % empty)
            self.tr_indentation ()
            if len (self._line) > 0 and \
                  self._line [0] == '#' and self._line [:2] != '##':
                self.tr_block_comment ()
                self._no_break = 1
            elif self._docprocess != 'none' and \
                 self._docstring and self._line[:3] in ('"""', "'''"):
                if self._docprocess == 'plain':
                    self.tr_docstring_plain()
                elif self._docprocess == 'struct':
                    self.tr_docstring_struct()
                else:
                    raise ValueError, 'Illegal value for doc process.'
            else:
                self._translation.append ('\\I{%d}{%d}' % \
                      (self._line_nr, self._indentation))
                self._docstuff.append ('\\I{%d}{%d}' % \
                      (self._line_nr, self._indentation))
                self.tr_code ()
                if not self._no_break and \
                      self._translation [-1] [-8:] == '\\colon $' and \
                      self._translation [0] [:3] != '\\E{':
                    self._translation.insert (0, '\\PB')
                    self._docstuff.insert (0, '\\PB')
                    self._no_break = 1
                else:
                    self._no_break = empty != 0
                if len (self._line) > 0:
                    if self._line [:1] != '#': raise error
                    if self._translation [-1] [:1] == '$':
                        self._translation.append ('\\quad ')
                        #self._docstuff.append ('\\quad ')
                    self.tr_comment ()
                else:
                    self._translation.append ('\n')
                    self._docstuff.append ('\n')
        except EOFError: pass
        return self._translation
    def translation (self):
        return self._docstuff
        #return self._translation

error = 'py2tex error'

class Re:
        def __init__(self, regex):
                self._regex = regex
        def match(self, string, pos = 0):
                m = self._regex.match(string, pos)
                result = -1
                if m:
                        result = m.end(0)
                return result
        def search(self, string, pos = 0):
                m = self._regex.search(string, pos)
                result = -1
                if m:
                        result = m.start(0)
                return result

class Regex:
        def compile(self, regex):
                return Re(re.compile(regex))

regex = Regex()

interpret_re = regex.compile ('%[ \t]*(La)?TeX[ \t]*$')
verbatim_re = regex.compile ('%[ \t]*ASCII[ \t]*$')
ordinary_re = regex.compile ('[^#$%&<>\\\\^_{|}~]*')
white_re = regex.compile ('[ \t]*')
function_re = regex.compile ('[ \t]*\\(')
comment_re = regex.compile ('(##|[^#])*')
escape_re = regex.compile ('([^\\\\]|\\\\.)*\\\\')
numeral_re = regex.compile (string.joinfields ((
      '0[xX][0-9A-Fa-f]+',
      '[0-9]+\\.?[eE][+-]?[0-9]+[jJLl]?',
      '[0-9]*\\.[0-9]+[eE][+-]?[0-9]+[jJLl]?',
      '[1-9][0-9]*[jJLl]?',
      '0[0-7]*'), '|'))

token_re = regex.compile (string.joinfields ((
      '[A-Za-z_][A-Za-z_0-9]*',
      "'('')?", '"("")?',
      '==?', '[<>!]=', '<>',
      '<<', '>>',
      '\\[]',
      '[*][*]',
      '[\\\\{}$&|^~%:*/+-]'), '|'))
quote_re= regex.compile( '("("")?)|' "('('')?)" ) 

TeX_code = {
      '\\': '$\\backslash$', '|': '$\\vert$',
      '<': '$<$', '>': '$>$',
      '{': '$\\{$', '}': '$\\}$' }
reserved = ('access', 'and', 'break', 'class', 'continue',
      'def', 'del', 'elif', 'else', 'except', 'exec',
      'finally', 'for', 'from', 'global', 'if',
      'import', 'in', 'is', 'is not', 'not', 'not in', 'or',
      'pass', 'print', 'raise', 'return', 'try', 'while')
single = ('else', 'finally', 'try', '-', '+')
reserved_operators = ('and', 'in', 'is', 'is not', 'not', 'not in', 'or', '**')
special_ctrl = { '\a': '\\a', '\b': '\\b', '\f': '\\f',
      '\n': '\\n', '\r': '\\r', '\t': '\\t', '\v': '\\v' }

def usc_protect (ident):
    ident = string.joinfields (string.splitfields (ident, '_'), '\\_')
    return ident

def ctrl_protect (str):
    result = ''
    for c in str:
        o = ord (c)
        if o < 32 or o >= 127:
            if special_ctrl.has_key (c):
                result = result + special_ctrl [c]
            else:
                result = '%s\\%03o' % (result, o)
        else:
            result = result + c
    return result

def tr_numeral (token):
    end = token[-1]                        # Preserve the type signifier (jJlL) if any.
    numeral = string.lower (token)
    if numeral [:2] == '0x':
        # #(0x1a, 0X2B)
        return '\\HEX{%s}' % string.upper (numeral [2:])
    if not (end in'jJlL'):                # Check if end is a signifier.
        end = ''
    else:
        numeral = numeral[:-1]                # Strip the signifier.
    pos = string.find (numeral, 'e')
    if pos >= 0:
        # #(12.4e-78, .3333E+0, .1e6, 2.e1, 0.e1, 1e4)
        return '\\EXP{%s}{%s}{%s}' % \
              (numeral [:pos], numeral [pos + 1:], end)
    if numeral [:1] == '0' and numeral != '0':
        # #(0377, 0378)
        return '\\OCT{%s}' % numeral [1:]
    # #(.333, 3.141592)# #(0, 1, 42)
    return '\\NUM{%s}{%s}' % (numeral, end) 

__author__ = 'Martin_Tiefengraber'

from numpy import *
from linecache import getline, clearcache
from sys import exit as sysexit


def GetLineEntry(filename, signal, delim, meth):
    # open text file
    fleParams = open("C:/Blick/src/CINDI3/" + filename, 'r')
    # localize wanted text line
    noline = None
    for i, l in enumerate(fleParams):
        if l.strip().partition(delim)[0] == signal:
            i += 1
            noline = False
            break
        else:
            noline = True
    # stop if no accordance was found
    if noline:
        print 'No entry "' + signal + '" in ' + filename + '!'
        sysexit(0)
    fleParams.close()
    # read wanted line
    entry = getline(filename, i).strip().partition(' -> ')[-1].split(';')
    clearcache()

    out = None
    # write wanted line to chosen Python format
    if meth == 'dictfloat':
        out = {}
        for i in entry:
            p = i.split(':')
            out[p[0].strip()] = asarray(map(float, p[1].split(',')))
    elif meth == 'dictint':
        out = {}
        for i in entry:
            p = i.split(':')
            out[p[0].strip()] = asarray(map(int, p[1].split(',')))
    elif meth == 'dictbool':
        out = {}
        for i in entry:
            p = i.split(':')
            if p[1].split(',')[0] == 'True':
                out[p[0].strip()] = asarray([True])
            elif p[1].split(',')[0] == 'False':
                out[p[0].strip()] = asarray([False])
            else:
                print 'unknown format in bool dict'
    elif meth == 'dictstring':
        out = {}
        for i in entry:
            p = i.split(':')
            out[p[0].strip()] = p[1].split(',')
    elif meth == 'float':
        if len(entry[0].split(',')) > 1:
            out = asarray(map(float, entry[0].split(',')))
        else:
            # out = array([entry[0]], dtype=float)
            out = array(entry[0], dtype=float)
            # out = float(entry[0])
    elif meth == 'int':
        if len(entry[0].split(',')) > 1:
            out = asarray(map(int, entry[0].split(',')))
        else:
            # out = array([entry[0]], dtype=int)
            out = array(entry[0], dtype=int)
    elif meth == 'string':
        if ',' in entry[0]:
            out = entry[0].split(',')
        else:
            out = entry[0]
    elif meth == 'stringarray':
        if ',' in entry[0]:
            out = asarray(entry[0].split(','))
        else:
            out = array(entry[0])
    elif meth == 'bool':
        if len(entry[0].split(',')) > 1:
            out = array([], dtype=bool)
            for e in entry[0].split(','):
                if e == 'True':
                    out = append(out, array(True, dtype=bool))
                elif e == 'False':
                    out = append(out, array(False, dtype=bool))
        else:
            if entry[0] == 'True':
                out = array(True, dtype=bool)
            elif entry[0] == 'False':
                out = array(False, dtype=bool)
    elif meth == 'dictbool':
        out = {}
        for i in entry:
            p = i.split(':')
            if p[1] == 'True':
                out[p[0].strip()] = array(True, dtype=bool)
            elif p[1] == 'False':
                out[p[0].strip()] = array(False, dtype=bool)

    return out
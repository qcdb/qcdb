import os
import sys


## {{{ http://code.activestate.com/recipes/52224/ (r1)
def search_file(filename, search_path):
    """Given an os.pathsep divided `search_path`, find first occurrence of
    `filename`. Returns full path to file if found or None if unfound.

    """
    file_found = False
    paths = search_path.split(os.pathsep)
    #paths = string.split(search_path, os.pathsep)
    for path in paths:
        if os.path.exists(os.path.join(path, filename)):
            file_found = True
            break
    if file_found:
        return os.path.abspath(os.path.join(path, filename))
    else:
        return None
## end of http://code.activestate.com/recipes/52224/ }}}


def all_casings(input_string):
    """Function to return a generator of all lettercase permutations
    of *input_string*.

    """
    if not input_string:
        yield ''
    else:
        first = input_string[:1]
        if first.lower() == first.upper():
            for sub_casing in all_casings(input_string[1:]):
                yield first + sub_casing
        else:
            for sub_casing in all_casings(input_string[1:]):
                yield first.lower() + sub_casing
                yield first.upper() + sub_casing


def import_ignorecase(module, lenv=None):
    """Function to import *module* in any possible lettercase
    permutation. Returns module object if available, None if not.
    `lenv` is list (not str) of addl sys.path members to try.

    """
    lenv = [] if lenv is None else lenv

    with add_path(lenv):
        modobj = None
        for per in list(all_casings(module)):
            try:
                modobj = __import__(per)
            except ImportError:
                pass
            else:
                break

    return modobj


class add_path():
    """https://stackoverflow.com/a/39855753"""

    def __init__(self, paths):
        # paths must be list
        self.paths = paths

    def __enter__(self):
        for pth in reversed(self.paths):
            sys.path.insert(0, pth)

    def __exit__(self, exc_type, exc_value, traceback):
        for pth in self.paths:
            sys.path.remove(pth)

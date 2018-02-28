import os

#def query_yes_no(question, default=True):
#    """Ask a yes/no question via raw_input() and return their answer.
#
#    *question* is a string that is presented to the user.
#    *default* is the presumed answer if the user just hits <Enter>.
#    It must be yes (the default), no or None (meaning
#    an answer is required of the user).
#
#    The return value is one of True or False.
#
#    """
#
#    yes = re.compile(r'^(y|yes|true|on|1)', re.IGNORECASE)
#    no = re.compile(r'^(n|no|false|off|0)', re.IGNORECASE)
#
#    if default == None:
#        prompt = " [y/n] "
#    elif default == True:
#        prompt = " [Y/n] "
#    elif default == False:
#        prompt = " [y/N] "
#    else:
#        raise ValueError("invalid default answer: '%s'" % default)
#
#    while True:
#        sys.stdout.write(question + prompt)
#        choice = raw_input().strip().lower()
#        if default is not None and choice == '':
#            return default
#        elif yes.match(choice):
#            return True
#        elif no.match(choice):
#            return False
#        else:
#            sys.stdout.write("    Please respond with 'yes' or 'no'.\n")


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


#def drop_duplicates(seq):
#    """Function that given an array *seq*, returns an array without any duplicate
#    entries. There is no guarantee of which duplicate entry is dropped.
#
#    """
#    #noDupes = []
#    #[noDupes.append(i) for i in seq if not noDupes.count(i)]
#    #return noDupes
#    noDupes = []
#    seq2 = sum(seq, [])
#    [noDupes.append(i) for i in seq2 if not noDupes.count(i)]
#    return noDupes
#
#
#def all_casings(input_string):
#    """Function to return a generator of all lettercase permutations
#    of *input_string*.
#
#    """
#    if not input_string:
#        yield ''
#    else:
#        first = input_string[:1]
#        if first.lower() == first.upper():
#            for sub_casing in all_casings(input_string[1:]):
#                yield first + sub_casing
#        else:
#            for sub_casing in all_casings(input_string[1:]):
#                yield first.lower() + sub_casing
#                yield first.upper() + sub_casing
#
#
#def getattr_ignorecase(module, attr):
#    """Function to extract attribute *attr* from *module* if *attr*
#    is available in any possible lettercase permutation. Returns
#    attribute if available, None if not.
#
#    """
#    array = None
#    for per in list(all_casings(attr)):
#        try:
#            getattr(module, per)
#        except AttributeError:
#            pass
#        else:
#            array = getattr(module, per)
#            break
#
#    return array
#
#
#def import_ignorecase(module):
#    """Function to import *module* in any possible lettercase
#    permutation. Returns module object if available, None if not.
#
#    """
#    modobj = None
#    for per in list(all_casings(module)):
#        try:
#            modobj = __import__(per)
#        except ImportError:
#            pass
#        else:
#            break
#
#    return modobj
#
#def findfile_ignorecase(fil, pre='', post=''):
#    """Function to locate a file *pre* + *fil* + *post* in any possible
#    lettercase permutation of *fil*. Returns *pre* + *fil* + *post* if
#    available, None if not.
#
#    """
#    afil = None
#    for per in list(all_casings(fil)):
#        if os.path.isfile(pre + per + post):
#            afil = pre + per + post
#            break
#        else:
#            pass
#
#    return afil


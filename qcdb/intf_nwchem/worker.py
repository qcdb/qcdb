import os
import uuid
import shutil
import subprocess


def nwchem_subprocess(nwchemrec):  # enginerec@i -> enginerec@io

    """
    Required Input Fields
    ---------------------
    command : list
        Command and arguments to execute. 
        Generally [xcfour]

    semioptional inputs b/c xmod

    Optional Input Fields
    ---------------------
    scratch_location : str, optional
        Override the default scratch location.
        Note that this IS ACTUAL (not PARENT) dir.
    executable_path
        Additional path (':'-separated if multiplie) to be prepended to `PATH` for subprocess execution

    Output Fields
    -------------
    stdout : str
        Main output file that gets written to stdout.
    output_grd : str, optional
        GRD file contents.
    output_fcmfinal : str, optional
    output_dipol : str, optional

    """

    try:
        nwchemrec['command']
        nwchemrec['nwchemnw']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            nwchemrec.keys())) from err

    current_directory = os.getcwd()

    # find environment by merging PSIPATH and PATH environment variables
    # * `path` kwarg gets precedence
    # * filter out None values as subprocess will fault on them
    lenv = {
        'HOME': os.environ.get('HOME'),
        'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                 ':' + os.environ.get('PATH')),# +
#                 ':' + qcdb.get_datadir() + '/basis'),
        'NWCHEM_OMP_NUM_CORES': os.environ.get('NWCHEM_OMP_NUM_CORES'),
#        'MKL_NUM_THREADS': os.environ.get('MKL_NUM_THREADS'),
#        'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    if 'executable_path' in nwchemrec:
        lenv['PATH'] = nwchemrec['executable_path'] + ':' + lenv['PATH']
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # set up unique scratch directory and move in
    if 'scratch_location' in nwchemrec:
        nwchem_tmpdir = nwchemrec['scratch_location']
    else:
        nwchem_tmpdir = os.environ['HOME'] + os.sep + 'nwchem_' + str(uuid.uuid4())
    if not os.path.exists(nwchem_tmpdir):
        os.mkdir(nwchem_tmpdir)
    os.chdir(nwchem_tmpdir)

    # write governing inputs
    with open('nwchem.nw', 'w') as handle:
        handle.write(nwchemrec['nwchemnw'])
#    with open('GENBAS', 'w') as handle:
#        handle.write(nwchemrec['genbas'])

    # call `xnwchem` program or subprogram
    try:
        spcall = subprocess.Popen(nwchemrec['command'], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as err:
        raise OSError('Command (`{}`) failed with PATH ({})'.format(
            ' '.join(nwchemrec['command']), lenv['PATH'])) from err

    # recover output data
    spcallstdout = ''
    while True:
        data = spcall.stdout.readline()
        data = data.decode('utf-8')
        if not data:
            break
#        print(data)
        spcallstdout += data
    nwchemrec['stdout'] = spcallstdout

#    for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
#        fullpath = nwchem_tmpdir + os.sep + fl
#        try:
#            with open(fullpath, 'r') as handle:
#                nwchemrec['output_' + fl.lower()] = handle.read()
#        except IOError:
#            pass

    # clean up files and remove scratch directory
    # NOTE used to keep scr arond if path in kwargs
    if 'scratch_messy' not in nwchemrec or nwchemrec['scratch_messy'] is False:
        os.chdir('..')
        try:
            shutil.rmtree(nwchem_tmpdir)
        except OSError as err:
            raise OSError('Unable to remove CFOUR temporary directory: {}'.
                          format(nwchem_tmpdir)) from err

    os.chdir(current_directory)

    return nwchemrec

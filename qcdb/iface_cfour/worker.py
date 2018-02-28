import os
import uuid
import shutil
import subprocess


def cfour_subprocess(cfourrec):  # cfourrec@i -> cfourrec@io

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
        cfourrec['command']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            cfourrec.keys())) from err

    current_directory = os.getcwd()

    # find environment by merging PSIPATH and PATH environment variables
    # * `path` kwarg gets precedence
    # * filter out None values as subprocess will fault on them
    lenv = {
        'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                 ':' + os.environ.get('PATH')),# +
#                 ':' + qcdb.get_datadir() + '/basis'),
#        'GENBAS_PATH': qcdb.get_datadir() + '/basis',
        'CFOUR_NUM_CORES': os.environ.get('CFOUR_NUM_CORES'),
        'MKL_NUM_THREADS': os.environ.get('MKL_NUM_THREADS'),
        'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    if 'executable_path' in cfourrec:
        lenv['PATH'] = cfourrec['executable_path'] + ':' + lenv['PATH']
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # set up unique scratch directory and move in
    if 'scratch_location' in cfourrec:
        cfour_tmpdir = cfourrec['scratch_location']
    else:
        cfour_tmpdir = os.environ['HOME'] + os.sep + 'cfour_' + str(uuid.uuid4())[:8]
    if not os.path.exists(cfour_tmpdir):
        os.mkdir(cfour_tmpdir)
    os.chdir(cfour_tmpdir)

    # write governing inputs
    with open('ZMAT', 'w') as handle:
        handle.write(cfourrec['zmat'])
    with open('GENBAS', 'w') as handle:
        handle.write(cfourrec['genbas'])

    # call `xcfour` program or subprogram
    try:
        spcall = subprocess.Popen(cfourrec['command'], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as err:
        raise OSError('Command (`{}`) failed with PATH ({})'.format(
            ' '.join(cfourrec['command']), lenv['PATH'])) from err

    # recover output data
    spcallstdout = ''
    while True:
        data = spcall.stdout.readline()
        data = data.decode('utf-8')
        if not data:
            break
#        print(data)
        spcallstdout += data
    cfourrec['stdout'] = spcallstdout

    for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
        fullpath = cfour_tmpdir + os.sep + fl
        try:
            with open(fullpath, 'r') as handle:
                cfourrec['output_' + fl.lower()] = handle.read()
        except IOError:
            pass

    # clean up files and remove scratch directory
    # NOTE used to keep scr arond if path in kwargs
    if 'scratch_messy' not in cfourrec or cfourrec['scratch_messy'] is False:
        #os.unlink(paramfileold)
        #os.unlink(paramfile)
        #os.unlink(geomfile)
        #if '-grad' in dftd3rec['command']:
        #    os.unlink(derivfile)

        os.chdir('..')
        try:
            shutil.rmtree(cfour_tmpdir)
        except OSError as err:
            raise OSError('Unable to remove CFOUR temporary directory: {}'.
                          format(cfour_tmpdir)) from err

    os.chdir(current_directory)

    return cfourrec

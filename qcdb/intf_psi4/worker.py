import os
import json
import uuid
import shutil
import subprocess


def psi4_subprocess(psi4rec):  # psi4rec@i -> psi4rec@io

    """
    Required Input Fields
    ---------------------
    command : list
        Command and arguments to execute.
        Generally [psi4 --json]
    json : json
        psi4 json input

    semioptional inputs b/c xmod

    Optional Input Fields
    ---------------------
#    scratch_location : str, optional
#        Override the default scratch location.
#        Note that this IS ACTUAL (not PARENT) dir.
    executable_path
        Additional path (':'-separated if multiplie) to be prepended to `PATH` for subprocess execution

    Output Fields
    -------------
#    stdout : str
#        Main output file that gets written to stdout.
#    output_grd : str, optional
#        GRD file contents.
#    output_fcmfinal : str, optional
#    output_dipol : str, optional

    """
    import bson

    try:
        psi4rec['command']
        psi4rec['json']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            psi4rec.keys())) from err

    current_directory = os.getcwd()
    jobuuid = str(uuid.uuid4())

    # set up unique scratch directory and move in
    if 'scratch_location' in psi4rec:
        tmpdir = psi4rec['scratch_location']
    else:
        tmpdir = os.environ['HOME'] + os.sep + 'psi4_' + jobuuid
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    os.chdir(tmpdir)

    # find environment by merging PSIPATH and PATH environment variables
    # * `path` kwarg gets precedence
    # * filter out None values as subprocess will fault on them
    lenv = {
        'HOME': os.environ.get('HOME'),
        'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                 ':' + os.environ.get('PATH')),# +
#                 ':' + qcdb.get_datadir() + '/basis'),
#        'GENBAS_PATH': qcdb.get_datadir() + '/basis',
#        'CFOUR_NUM_CORES': os.environ.get('CFOUR_NUM_CORES'),
#        'MKL_NUM_THREADS': os.environ.get('MKL_NUM_THREADS'),
#        'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
        'PSI_SCRATCH': tmpdir,
        'PYTHONPATH': os.environ.get('PYTHONPATH'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    if 'executable_path' in psi4rec:
        lenv['PATH'] = psi4rec['executable_path'] + ':' + lenv['PATH']
    lenv = {k: v for k, v in lenv.items() if v is not None}


    # write governing inputs
    inputjson = jobuuid + '.json'
    with open(inputjson, 'w') as handle:
        json.dump(psi4rec['json'], handle)
    #with open(inputjson, 'wb') as handle:
    #    handle.write(bson.dumps(psi4rec['json']))
    psi4rec['command'].append(inputjson)

    for fl in psi4rec['json'].keys():
        if fl.startswith('infile_'):
            with open(fl[7:], 'w') as handle:
                handle.write(psi4rec['json'][fl])

    # call `psi4` program
    try:
        spcall = subprocess.Popen(psi4rec['command'], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as err:
        raise OSError('Command (`{}`) failed with PATH ({})'.format(
            ' '.join(psi4rec['command']), lenv['PATH'])) from err

    # recover output data
    spcallstdout = ''
    while True:
        data = spcall.stdout.readline()
        data = data.decode('utf-8')
        if not data:
            break
        #print(data)
        spcallstdout += data
    psi4rec['stdout'] = spcallstdout

    #with open(inputjson, 'rb') as handle:
    #    psi4rec['json'] = bson.loads(handle.read())
    with open(inputjson, 'r') as handle:
        psi4rec['json'] = json.load(handle)

    for fl in ['grid_esp.dat']:
        fullpath = tmpdir + os.sep + fl
        try:
            with open(fullpath, 'r') as handle:
                psi4rec['json']['outfile_' + fl.lower()] = handle.read()
                print(handle.read())
        except IOError:
            pass

    # clean up files and remove scratch directory
    # NOTE used to keep scr arond if path in kwargs
    if 'scratch_messy' not in psi4rec or psi4rec['scratch_messy'] is False:
        os.chdir('..')
        try:
            shutil.rmtree(tmpdir)
        except OSError as err:
            raise OSError('Unable to remove Psi4 temporary directory: {}'.
                          format(tmpdir)) from err

    os.chdir(current_directory)


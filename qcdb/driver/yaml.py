from ..molecule import Molecule

def yaml_run(yamlin):

    import yaml
    cmd = yaml.load(yamlin)
    print('CMD', cmd)
    kwargs = cmd.get('kwargs', {})

    _, jrec = cmd['driver'](cmd['method'],
                            options=cmd['options'],
                            molecule=Molecule(cmd['molecule']),
                            return_wfn=True,
                            **kwargs)
                            #**cmd['kwargs'])
    return jrec


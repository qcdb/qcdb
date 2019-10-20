from ..molecule import Molecule


def yaml_run(yamlin):

    import yaml
    cmd = yaml.load(yamlin, Loader=yaml.FullLoader)

    kwargs = cmd.get('kwargs', {})

    _, jrec = cmd['driver'](cmd['method'],
                            options=cmd['options'],
                            molecule=Molecule(cmd['molecule']),
                            return_wfn=True,
                            **kwargs)
    return jrec

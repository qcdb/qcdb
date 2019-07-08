
def gamess_list():
        """Return an array of gamess methods with energies. Appended
    to procedures['energy'].
        """
        val = []
        val.append('gamess')
        val.append('gms-gamess')
        val.append('gms-makefp')
        val.append('gms-efp')
        val.append('gms-scf')
        val.append('gms-hf')
        val.append('gms-mp2')
        val.append('gms-ccsd')
        val.append('gms-ccsd(t)')
        val.append('gms-ccsd(tq)')
        val.append('gms-fci')
        val.append('gms-b3lyp')

        return val


def gamess_gradient_list():
        """Return an array of gamess methods with energies. Appended
    to procedures['gradient'].
        """
        val = []
        return val

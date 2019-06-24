
def gamess_list():
        """Return an array of Cfour methods with energies. Appended
    to procedures['energy'].
        """
        val = []
        val.append('gamess')
        val.append('gms-scf')
        val.append('gms-hf')
        val.append('gms-mp2')
        val.append('gms-ccsd')
        val.append('gms-ccsd(t)')
        val.append('gms-ccsd(tq)')
        val.append('gms-fci')
#        val.append('gms-dft')
#        val.append('gms-efp')

        return val


def gamess_gradient_list():
        """Return an array of Cfour methods with energies. Appended
    to procedures['gradient'].
        """
        val = []
        val.append('gms-hf')
        val.append('gms-scf')
        return val

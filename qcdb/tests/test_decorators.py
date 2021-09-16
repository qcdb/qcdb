import numpy as np

import qcdb
from qcdb.keywords.decorators import def_mol

mol1 = qcdb.Molecule(
    """
He
"""
)

mol3 = qcdb.Molecule(
    """
He
He 1 1.5
He 1 1.5 2 90.0
"""
)

mol4 = qcdb.Molecule(
    """
He 0 0 0
He 1 1 1
He 2 2 2
He 3 3 3
"""
)

mol0 = qcdb.Molecule("""Xe""")


@def_mol()
def fenergy(name, **kwargs):
    print(kwargs["molecule"])
    return kwargs["molecule"].natom()


@def_mol()
def fgrad(name, **kwargs):
    limol = kwargs["molecule"].natom()
    assert fenergy(name, molecule=kwargs["molecule"]) == limol
    numol = qcdb.Molecule(geom=np.arange(limol * 6).reshape((-1, 3)), elez=[2 for _ in range(limol * 2)])
    numol.update_geometry()
    assert fenergy(name, molecule=numol) == limol * 2
    return numol.natom()


# can't test en masse b/c P::e not fresh
# def test_1():
#    # from default H2 in P::e
#    assert fenergy('mbpt') == 2


def test_2():
    qcdb.activate(mol3)
    assert fenergy("mbpt") == 3


def test_3():
    assert fenergy("mbpt", molecule=mol4) == 4


def test_4():
    qcdb.set_molecule(
        """
He
He 1 1.5
He 1 1.5 2 90.0
"""
    )
    assert fenergy("mbpt") == 3


def test_5():
    assert fgrad("uuio", molecule=mol4) == 8


def test_6():
    qcdb.activate(mol1)
    assert fgrad("uuio") == 2


def hide_test_24():
    subjects = RottenOptions()
    subjects.add("qcdb", RottenOption(keyword="scf_e_conv", default=5, validator=parsers.parse_convergence))

    @register_opts(subjects)
    def energy(count=1, **kwargs):
        if count > 3:
            assert subjects.scroll["QCDB"]["SCF_E_CONV"].value == 1.0e-4
            return
        else:
            count += 1
            # subjects.require('qcdb', 'scf_c_conv', count, accession=kwargs['accession'])
            subjects.require("qcdb", "scf_e_conv", count, accession=kwargs["accession"])
            proc(count)

    @register_opts(subjects)
    def proc(count, **kwargs):
        energy(count)

    assert subjects.scroll["QCDB"]["SCF_E_CONV"].value == 1.0e-5
    energy()
    assert subjects.scroll["QCDB"]["SCF_E_CONV"].value == 1.0e-5

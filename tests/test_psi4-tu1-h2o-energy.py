import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *


@using_psi4
def test_tu1a():
    """tu1-h2o-energy/input.dat
    global testing

    """
#    import qcdb
#    print('qcdb', qcdb.__file__)
#    print(dir(qcdb))
    h2o = qcdb.set_molecule("""
  O 
  H 1 0.96
  H 1 0.96 2 104.5
""")
    print(h2o)
    print(qcdb.get_active_options().print_changed())

    qcdb.set_options({'basis': 'cc-pVDZ',
                      'memory': '600 mb'})

    qcdb.energy('scf')
    print(qcdb.print_variables())

    #assert compare_values(-76.0266327341067125, get_variable('SCF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert compare_values(-76.0266327341067125, qcdb.get_variable('SCF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
#    assert False


@using_psi4
def test_tu1b():
    """tu1-h2o-energy/input.dat
    local testing

    """

#    memory 600 mb

    h2o = qcdb.Molecule("""
  O 
  H 1 0.96
  H 1 0.96 2 104.5
""")

    E, jrec = qcdb.energy('scf/cc-pVDZ', molecule=h2o, return_wfn=True)
    print(qcdb.print_variables(jrec['qcvars']))

    ans1 = -76.0266327341067125
    assert compare_values(ans1, jrec['qcvars']['SCF TOTAL ENERGY'].data, 6, sys._getframe().f_code.co_name)
    # ok in direct mode
    #assert compare_values(ans1, jrec['wfn'].get_variable('SCF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)


@using_psi4
def test_tu2():
    """tu2-ch2-energy/input.dat
    #! Sample UHF/6-31G** CH2 computation

    """
    ans2 = -38.9253416082900827


    ch2 = qcdb.set_molecule("""
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93
""")

    qcdb.set_options({'basis': '6-31G**',
                      'reference': ' uhf',
                      'puream': 'cart',
                      #'psi_scf_type': 'pk'})
                      'scf_type': 'pk'})

    E, jrec = qcdb.energy ('scf', return_wfn=True)
    print(qcdb.print_variables())

    assert compare_values(ans2, qcdb.get_variable('scf total energy'), 6, sys._getframe().f_code.co_name)


#@using_psi4
#def test_tu2_sowreap():
#    """tu2-ch2-energy/input.dat
#    #! Sample UHF/6-31G** CH2 computation
#
#    """
#    ans2 = -38.9253416082900827
#
#
#    ch2 = qcdb.set_molecule("""
#  0 3
#  C
#  H 1 R
#  H 1 R 2 A
#
#  R = 1.075
#  A = 133.93
#""")
#
#    qcdb.set_options({'basis': '6-31G**',
#                      'reference': ' uhf',
#                      'puream': 'cart',
#                      #'psi_scf_type': 'pk'})
#                      'scf_type': 'pk'})
#
#    E, jrec = qcdb.energy ('scf', return_wfn=True, probe=True)
#    print(qcdb.print_variables())
#
#    assert compare_values(ans2, qcdb.get_variable('scf total energy'), 6, sys._getframe().f_code.co_name)

@using_psi4
def test_tu2_yaml():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: scf
options:
  memory: 1gb
  basis: '6-31g**'
  reference: uhf
  puream: cart
  scf_type: pk
"""

    import yaml
    asdf = yaml.load(yamlin)

    ene = asdf['driver'](asdf['method'],
                         options=asdf['options'],
                         molecule=qcdb.Molecule(asdf['molecule']))

    assert compare_values(-38.9253416082900827, ene, 6, 'calc from yaml str')





if __name__ == '__main__':
    test_tu1a()
    #test_tu1b()
    #test_tu2()

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import qcdb

_ref_h2o_pk_rhf  = -76.02696997325441
_ref_ch2_pk_uhf  = -38.925088643363665
_ref_ch2_pk_rohf = -38.91973113928147

@using_gamess
def test_tu1_rhf_a():
    """tu1-h2o-energy/input.dat
    global testing

    """
    h2o = qcdb.set_molecule("""
  O 
  H 1 1.8
  H 1 1.8 2 104.5
  units au
""")
    print(h2o)
    print(qcdb.get_active_options().print_changed())

    qcdb.set_options({'basis': 'cc-pVDZ',
                      'memory': '600 mb'})

    qcdb.energy('gms-hf')
    print(qcdb.print_variables())

    assert compare_values(_ref_h2o_pk_rhf, qcdb.variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)


@using_gamess
def test_tu1_rhf_b():
    """tu1-h2o-energy/input.dat
    local testing

    """

#    memory 600 mb

    h2o = qcdb.Molecule("""
  O 
  H 1 1.8
  H 1 1.8 2 104.5
  units au
""")

    E, jrec = qcdb.energy('gms-hf/cc-pVDZ', molecule=h2o, return_wfn=True)
    print(qcdb.print_variables(jrec['qcvars']))

    assert compare_values(_ref_h2o_pk_rhf, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, sys._getframe().f_code.co_name)


@using_gamess
def test_tu2_uhf():
    """tu2-ch2-energy/input.dat
    #! Sample UHF/6-31G** CH2 computation

    """


    ch2 = qcdb.set_molecule("""
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 2.05
  A = 133.93
units au
""")

    qcdb.set_options({'basis': '6-31G**',
                      'reference': ' uhf',
                      'puream': 'cart',
                      #'psi_scf_type': 'pk'})
                      'scf_type': 'pk'})

    E, jrec = qcdb.energy ('gms-hf', return_wfn=True)
    print(qcdb.print_variables())

    assert compare_values(_ref_ch2_pk_uhf, qcdb.variable('hf total energy'), 6, sys._getframe().f_code.co_name)


@using_gamess
def test_tu2_rohf():
    """tu2-ch2-energy/input.dat
    #! Sample ROHF/6-31G** CH2 computation

    """

    ch2 = qcdb.set_molecule("""
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 2.05
  A = 133.93
units au
""")

    qcdb.set_options({'basis': '6-31G**',
                      'reference': ' rohf',
                      'puream': 'cart',
                      #'psi_scf_type': 'pk'})
                      'scf_type': 'pk'})

    E, jrec = qcdb.energy ('gms-hf', return_wfn=True)
    print(qcdb.print_variables())

    assert compare_values(_ref_ch2_pk_rohf, qcdb.variable('hf total energy'), 6, sys._getframe().f_code.co_name)
    assert compare_values(_ref_ch2_pk_rohf, qcdb.variable('current energy'), 6, sys._getframe().f_code.co_name)
    assert compare_values(_ref_ch2_pk_rohf, E, 6, sys._getframe().f_code.co_name)


#@using_gamess
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
#    E, jrec = qcdb.energy ('gms-scf', return_wfn=True, probe=True)
#    print(qcdb.print_variables())
#
#    assert compare_values(ans2, qcdb.variable('scf total energy'), 6, sys._getframe().f_code.co_name)

@using_gamess
def test_tu2_uhf_yaml():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: gms-hf
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

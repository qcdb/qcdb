import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from utils import *
from addons import *


sys2_pk_uhf =  -38.925334624580
sys2_pk_rohf = -38.919988755499

@using_nwchem
def test_tu1a():
    """tu1-h2o-energy/input.dat
    global testing

    """
    h2o = qcdb.set_molecule("""
  O 
  H 1 0.96
  H 1 0.96 2 104.5
""")
    print(h2o)
    print(qcdb.get_active_options().print_changed())

    qcdb.set_options({'basis': 'cc-pVDZ',
                      'memory': '600 mb'})

    qcdb.energy('nwc-hf')
    print(qcdb.print_variables())

    pk_ref = -76.0266537
    assert compare_values(pk_ref, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)


@using_nwchem
def test_tu1b():
    """tu1-h2o-energy/input.dat
    local testing

    """

    h2o = qcdb.Molecule("""
  O 
  H 1 0.96
  H 1 0.96 2 104.5
""")

    E, jrec = qcdb.energy('nwc-hf/cc-pVDZ', molecule=h2o, return_wfn=True)
    print(qcdb.print_variables(jrec['qcvars']))

    #ans1 = -76.0266327341067125
    pk_ref = -76.0266537
    assert compare_values(pk_ref, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, sys._getframe().f_code.co_name)
    assert compare_values(pk_ref, jrec['qcvars']['SCF TOTAL ENERGY'].data, 6, sys._getframe().f_code.co_name)


@using_nwchem
def test_tu2_rohf_api():
    """tu2-ch2-energy/input.dat
    #! Sample UHF/6-31G** CH2 computation

    """
    ans2 = -38.9199888 #uhf-38.9253416082900827

    # ans2 = -38.91999637845108  # psi pk rohf
    # ans2 = -38.92534160829007  # psi pk uhf

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

    E, jrec = qcdb.energy ('nwc-scf', return_wfn=True)
    print(qcdb.print_variables())

    # nwc not well converged. tighten after we can dial in by option
    assert compare_values(ans2, qcdb.get_variable('scf total energy'), 4, sys._getframe().f_code.co_name)


@using_nwchem
def test_tu2_rohf_yaml():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: nwc-hf
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

    # UHF assert compare_values(-38.9253416082900827, ene, 6, 'calc from yaml str')
    assert compare_values(-38.9199888, ene, 4, 'calc from yaml str')  # ROHF





if __name__ == '__main__':
    test_tu1a()
    #test_tu1b()
    #test_tu2()

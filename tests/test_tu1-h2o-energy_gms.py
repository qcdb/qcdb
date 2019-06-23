import os
import pprint as pp
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

    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)

@using_gamess
def test_tu1_mp2():
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

    qcdb.energy('gms-mp2')
    print(qcdb.print_variables())

    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0

@using_gamess
def test_tu1_dft():
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

    qcdb.energy('gms-b3lyp')
    print(qcdb.print_variables())

    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0

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

    assert compare_values(_ref_ch2_pk_uhf, qcdb.get_variable('hf total energy'), 6, sys._getframe().f_code.co_name)


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

    assert compare_values(_ref_ch2_pk_rohf, qcdb.get_variable('hf total energy'), 6, sys._getframe().f_code.co_name)
    assert compare_values(_ref_ch2_pk_rohf, qcdb.get_variable('current energy'), 6, sys._getframe().f_code.co_name)
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
#    assert compare_values(ans2, qcdb.get_variable('scf total energy'), 6, sys._getframe().f_code.co_name)

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



@using_gamess
def test_makefp_h2o_h2o():
    """makefp_h2o/input.dat
    global testing

    """
    h2o = qcdb.set_molecule("""
  O 
  H 1 1.8
  H 1 1.8 2 104.5
  units au
""")

    qcdb.set_options({'basis': 'cc-pVDZ',
                      'memory': '600 mb',})
                      #'gamess_contrl__runtyp': 'makefp'})

    E, jobrec = qcdb.energy('gms-makefp', return_wfn=True) #gamess')

    efppot = jobrec['extras']['fragment_potential']
    mysteryvalue = len(jobrec['molecule']['symbols'])


    h2o_h2o = qcdb.set_molecule("""
 O        -2.47693896     1.68753624     0.00297498
 H        -2.47693896     2.44253635    -0.59192497
 H        -2.47693896     0.93253607    -0.59192497
 O         0.99148941     0.37922010     0.00256592
 H         0.99148941     1.13422012    -0.59233409
 H         0.99148941    -0.37577990    -0.59233409
 O        -3.01032066    -1.88144445    -0.01070227
 H        -3.01032066    -1.12644446    -0.60560226
 H        -3.01032066    -2.63644457    -0.60560226
 O         1.24186850    -2.66617513    -0.00931387
 H         1.24186850    -1.91117513    -0.60421389
 H         1.24186850    -3.42117500    -0.60421389
 O         3.82095551    -1.43027246    -0.00164440
 H         3.82095551    -0.67527246    -0.59654438
 H         3.82095551    -2.18527246    -0.59654438
  units au
""")

#    print(h2o_h2o)

    print(qcdb.get_active_options().print_changed())

    qcdb.set_options({'basis': 'cc-pVDZ',
                      'memory': '600 mb'})
                      #'gamess_contrl__coord': 'unique'})

#    E, jobrec = qcdb.energy('gms-efp', return_wfn=True, efpfrag=efppot) #gamess')

    E, jobrec = qcdb.energy('gms-efp', return_wfn=True, efpfrag=efppot, number_of_atoms_in_frag_x=mysteryvalue) #gamess')

    print(qcdb.print_variables())

    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0



@using_gamess
def test_makefp_h2o_nh3():
    """makefp_h2o/input.dat
    global testing

    """
    h2o = qcdb.set_molecule("""
 O              6.3036714007  -2.4307101077   0.4363340550
 H              6.7562107152  -1.6825574097   0.7808716709
 H              6.8187921909  -3.1820001373   0.6675115850
  units au
""")

    qcdb.set_options({'basis': '6-31G*',
                      'memory': '600 mb',
                      'makefp__frag' : 'water'})

    E, jobrec = qcdb.energy('gms-makefp', return_wfn=True)

    efppot_h2o = jobrec['extras']['fragment_potential']
    print(efppot_h2o)

    nh3 = qcdb.set_molecule("""
N        -0.27418831    -0.27579930     2.13535595
H        -0.27418831     0.79420060     2.13535595
H        -0.27418831    -0.63246512     3.14416194
H         0.59946275    -0.63246512     1.63095403
  units au
""")

    qcdb.set_options({'basis': '6-31G*',
                      'memory': '600 mb',})

    E, jobrec = qcdb.energy('gms-makefp', return_wfn=True) 

    efppot_nh3 = jobrec['extras']['fragment_potential']
    print(efppot_nh3)


    h2o_nh3 = qcdb.set_molecule("""
O        -0.49406338    -1.51857603     0.12618840
H        -0.17619392    -0.87140638     0.76184309
H        -0.08951682    -2.33959365     0.41978461
N         0.55317599    -0.10685764     2.13535643
H         0.55317599     0.96314228     2.13535643
H         0.55317599    -0.46352431     3.14416218
H         1.42682731    -0.46352422     1.63095355
  units au
""")

    qcdb.set_options({'basis': '6-31G*',
                      'memory': '600 mb'})

    E, jobrec = qcdb.energy('gms-efp', return_wfn=True, efpfrag1=efppot_h2o, efpfrag2=efppot_nh3)

    pp.pprint(qcdb.print_variables())

    print(jobrec)

    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0



ref  = 'Cnv'

@using_gamess
def test_tu1_symmetry():
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

    print(qcdb.get_variable('THE POINT GROUP OF THE MOLECULE'))
#    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('THE POINT GROUP OF THE MOLECULE'), 6, sys._getframe().f_code.co_name)
    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0


if __name__ == '__main__':
    test_tu1a()
    #test_tu1b()
    #test_tu2()

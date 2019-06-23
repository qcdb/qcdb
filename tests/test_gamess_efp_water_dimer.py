import os
import pprint as pp
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import qcdb

#_ref_h2o_pk_rhf  = -76.02696997325441
#_ref_ch2_pk_uhf  = -38.925088643363665
#_ref_ch2_pk_rohf = -38.91973113928147

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

#    assert compare_values(_ref_h2o_pk_rhf, qcdb.get_variable('HF TOTAL ENERGY'), 6, sys._getframe().f_code.co_name)
    assert 0


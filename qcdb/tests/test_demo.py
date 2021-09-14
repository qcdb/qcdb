import pprint
import qcdb
import psi4
pp = pprint.PrettyPrinter(width=120)

""" Helper Functions """


def gen_rvals(center, stepsize , npoints):
  if npoints % 2 == 0:
    return []
  center_ind = npoints//2
  rvals = [center + stepsize*(ind - center_ind) for ind in range(npoints)]
  return rvals


""" QCDB Setup """

hf = qcdb.set_molecule("""
        B
        H 1 R
        units angstrom
        """)

local_options = {"memory": 35, "nnodes": 1, "ncores": 1}
qcdb.set_options({
        'e_convergence': 1e-11,
        'scf__d_convergence': 1e-9,
        'nwchem_ccsd__maxiter': 100,
        'psi4_mp2_type': 'conv',
        'psi4_scf_type': 'direct',
        'psi4_df_scf_guess': 'false',
                    })


""" Geometry Scan """

Re = 1.229980860371199 # self consistent with fci optimization
dR = 0.005
npoints = 5

R_arr = gen_rvals(Re, dR, npoints)

# energy of base method
E_base = [0.0] * npoints

# various corrections to base energy
dE_basis = [0.0] * npoints
dE_dboc = [0.0] * npoints
dE_x2c = [0.0] * npoints
dE_fci = [0.0] * npoints
dE_ccsdtq = [0.0] * npoints

for i, R in enumerate(R_arr):
  print(f'<<<\n<<< POINT {i}\n<<<')

  hf.R = R
  print(f'~~~ Calculation at Re={R} Angstrom ({i+1}/{npoints}) ~~~')

  # read from file in case of restart
#  if i < 1:
#  if True:
#  if i > 0:
  if False:
#  if i < 4:
      with open(f'dist_{i}', 'r+') as f:
          lines = f.readlines()
          lines = [float(line.split()[3]) for line in lines]
          E_base[i] = lines[0]
          dE_basis[i] = lines[1]
          dE_dboc[i] = lines[2]
          dE_x2c[i] = lines[3]
          dE_fci[i] = lines[4]
          dE_ccsdtq[i] = lines[5]
          print(lines)
      continue

  # ccsdtq correction: (CCSDTQ - CCSD(T)) / cc-pVDZ
  qcdb.set_options({'cfour_dropmo': [1],})
  _, jrec = qcdb.energy('c4-ccsd(t)/cc-pVTZ', return_wfn=True, local_options=local_options)
  E_ccsdpt = float(jrec['qcvars']['CCSD(T) TOTAL ENERGY'].data)
  _, jrec = qcdb.energy('c4-ccsdtq/cc-pVTZ', return_wfn=True, local_options=local_options)
  E_ccsdtq = float(jrec['qcvars']['CCSDTQ TOTAL ENERGY'].data)
  qcdb.set_options({'cfour_dropmo': None})
  dE_ccsdtq[i] = E_ccsdtq - E_ccsdpt
  print(f'~~~ CCSDTQ Correction={dE_ccsdtq[-1]} Har. ({i+1}/{npoints}) ~~~')

  # base calculation: CCSD(T) / cc-pCV[Q5]Z
#  qcdb.set_options({'memory': '10 gb'})
  #E, jrec = qcdb.energy('nwc-ccsd(t)/cc-pCV[T,Q]Z', return_wfn=True)
  E, jrec = qcdb.energy('nwc-ccsd(t)/cc-pCVTZ', return_wfn=True)
#  qcdb.set_options({'memory': '55 gb'})
  E_base[i] = E
  print(f'~~~ Base Energy={E} Har. ({i+1}/{npoints}) ~~~')

  # basis set correction: MP2 / (aug-cc-pCV[56]Z) - cc-pCV[Q5]Z)
  E_small, _ = qcdb.energy('p4-mp2/cc-pCV[T,Q]Z', return_wfn=True)
  E_large, _ = qcdb.energy('p4-mp2/aug-cc-pCV[T,Q]Z', return_wfn=True)
  dE_basis[i] = E_large - E_small
  print(f'~~~ Basis Correction={dE_basis[-1]} Har. ({i+1}/{npoints}) ~~~')

  # relativistic correction: (X2C-CCSD(T) - CCSD(T)) / cc-pCVTZ-DK
  qcdb.set_options({'psi4_relativistic': 'x2c'})
  E_x2c_on, jrec = qcdb.energy('p4-ccsd(t)/aug-cc-pCVTZ-DK', return_wfn=True)
  qcdb.set_options({'psi4_relativistic': 'no'})
  E_x2c_off, jrec = qcdb.energy('p4-ccsd(t)/aug-cc-pCVTZ-DK', return_wfn=True)
  dE_x2c[i] = E_x2c_on-E_x2c_off
  print(f'~~~ Relativistic Correction={dE_x2c[-1]} Har. ({i+1}/{npoints}) ~~~')

  # fci correction: (FCI - CCSD(T)) / cc-pVDZ
  E_cc, _ = qcdb.energy('gms-ccsd(t)/cc-pVDZ', return_wfn=True, local_options=local_options)
  E_fci, _ = qcdb.energy('gms-fci/cc-pVDZ', return_wfn=True, local_options=local_options)
  dE_fci[i] = E_fci - E_cc
  print(f'~~~ FCI Correction={dE_fci[-1]} Har. ({i+1}/{npoints}) ~~~')

  with open(f'dist_{i}', 'a+') as f:
    f.write(f'Re {R_arr[i]} ,E_base     {E_base[i]} ')
    f.write('\n')
    f.write(f'Re {R_arr[i]} ,dE_basis   {dE_basis[i]} ')
    f.write('\n')
    f.write(f'Re {R_arr[i]} ,dE_dboc    {dE_dboc[i]} ')
    f.write('\n')
    f.write(f'Re {R_arr[i]} ,dE_x2c     {dE_x2c[i]} ')
    f.write('\n')
    f.write(f'Re {R_arr[i]} ,dE_fci     {dE_fci[i]} ')
    f.write('\n')
    f.write(f'Re {R_arr[i]} ,dE_ccsdtq  {dE_ccsdtq[i]} ')
    f.write('\n')

  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
  print(f'Re {R_arr[i]} ,E_base     {E_base[i]}')
  print(f'Re {R_arr[i]} ,dE_basis   {dE_basis[i]}')
  print(f'Re {R_arr[i]} ,dE_x2c     {dE_x2c[i]}')
  print(f'Re {R_arr[i]} ,dE_fci     {dE_fci[i]}')
  print(f'Re {R_arr[i]} ,dE_ccsdtq  {dE_ccsdtq[i]}')
  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

# Energy with a single correction
E_basis = []
E_dboc = [] 
E_x2c = [] 
E_fci = [] 
E_ccsdtq = []
# Total energies, using all corrections
E_tot_fci = []
E_tot_ccsdtq = []
for i in range(npoints):

    E_basis.append(E_base[i] + dE_basis[i])
    E_dboc.append(E_base[i] + dE_dboc[i])
    E_x2c.append(E_base[i] + dE_x2c[i])
    E_fci.append(E_base[i] + dE_fci[i])
    E_ccsdtq.append(E_base[i] + dE_ccsdtq[i])

    E_tot_fci.append(E_base[i] + dE_basis[i] + dE_dboc[i] + dE_x2c[i] + dE_fci[i])
    E_tot_ccsdtq.append(E_base[i] + dE_basis[i] + dE_dboc[i] + dE_x2c[i] + dE_ccsdtq[i])

for i in range(len(R_arr)):
    print('Re', R_arr[i], )
    print('   E_base', E_base[i])
    print('   dE_basis', dE_basis[i])
    print('   dE_x2c', dE_x2c[i])
    print('   dE_fci', dE_fci[i])
    print('   dE_ccsdtq', dE_ccsdtq[i])
    print('   E_tot_fci', E_tot_fci[i])
    print('   E_tot_ccsdtq', E_tot_ccsdtq[i])

""" Spectroscopic Constants """

# need a Psi4 molecule to give atomic masses to the diatomic module
psi4.geometry("""
        B
        H 1 1000
        units au
        """)

print('Total correction w/ FCI')
phys_consts_tot_fci = psi4.diatomic.anharmonicity(R_arr, E_tot_fci)
pp.pprint(phys_consts_tot_fci)
print('Total correction w/ CCSDTQ')
phys_consts_tot_ccsdtq = psi4.diatomic.anharmonicity(R_arr, E_tot_ccsdtq)
pp.pprint(phys_consts_tot_ccsdtq)

print('Base Energy')
phys_consts_base = psi4.diatomic.anharmonicity(R_arr, E_base)
pp.pprint(phys_consts_base)

print('Base Energy w/ only basis correction')
phys_consts_basis = psi4.diatomic.anharmonicity(R_arr, E_basis)
pp.pprint(phys_consts_basis)
print('Base Energy w/ only DBOC')
phys_consts_dboc = psi4.diatomic.anharmonicity(R_arr, E_dboc)
pp.pprint(phys_consts_dboc)
print('Base Energy w/ only rel. correction')
phys_consts_x2c = psi4.diatomic.anharmonicity(R_arr, E_x2c)
pp.pprint(phys_consts_x2c)
print('Base Energy w/ only FCI correction')
phys_consts_fci = psi4.diatomic.anharmonicity(R_arr, E_fci)
pp.pprint(phys_consts_fci)
print('Base Energy w/ only CCSDTQ correction')
phys_consts_ccsdtq = psi4.diatomic.anharmonicity(R_arr, E_ccsdtq)
pp.pprint(phys_consts_ccsdtq)



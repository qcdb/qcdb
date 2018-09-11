#! single-point UHF-MP2/cc-pVDZ  on NH2 
# ROHF-MP2 is not available in NWChem 
print('        <<< Literal nwchem.nw to NWChem >>>')

print('''memory 300 mb

molecule {
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
}

nwchem {
basis spherical
* library cc-pVDZ
end

scf
uhf
nopen 1
thresh 1.0e-8
maxiter 80
end

task mp2 energy
}

energy('nwchem')


clean()
clean_variables()
nwchem {}
''')

nh2= qcdb.set_molecule('''
        0 2
        N
        H 1 R
        H 1 1.008 2 105.0''')

print(nh2)
qcdb.set_options({
    'basis': 'cc-pvdz',
    'memory': '300 mb',
    'nwchem_scf': 'UHF',
    'nwchem_scf_nopen': 1,
    'nwchem_scf_thresh': 1.0e-8,
    'nwchem_scf_maxiter': 80,
    'nwchem_task_mp2': 'energy'
    })
print(qcdb.get-active_options().print_changed())

def check_rohf_mp2(return_value, is_df):
    if is_df:
        ref=
        mp2=


print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>'

banner('UHF-MP2 energy calculation')

molecule {
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
}

set {
basis cc-pVDZ
nwchem_scf uhf
nwchem_scf_nopen 1
nwchem_scf_maxiter 80
nwchem_scf_thresh 1.0e-8
}

energy('nwchem-mp2')


clean()
clean_variables()
nwchem {}
''')
print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>

banner('UHF-MP2 energy calculation')

molecule {
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
}

set {
basis cc-pvdz
reference uhf
maxiter 80 
nwchem_scf_nopen 1
nwchem_scf_thresh 1.0e-8
}

energy('nwchem-mp2')

''')

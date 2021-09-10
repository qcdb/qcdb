# <img src="https://github.com/qcdb/qcdb/blob/master/media/qcdb_newlogo.png" height=150>

| **Status** | [![GHA build](https://img.shields.io/github/workflow/status/qcdb/qcdb/CI?logo=github)](https://github.com/qcdb/qcdb/actions/workflows/ci.yml) [![Codecov coverage](https://img.shields.io/codecov/c/github/qcdb/qcdb.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/qcdb/qcdb) [![LGTM analysis](https://img.shields.io/lgtm/grade/python/g/qcdb/qcdb.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/qcdb/qcdb/context:python) |
| :------ | :------- |
| **Communication** | [![docs latest](https://img.shields.io/badge/docs-latest-5077AB.svg?logo=read%20the%20docs)](https://qcdb.github.io/qcdb/) [![dev chat on slack](https://img.shields.io/badge/dev_chat-on_slack-808493.svg?logo=slack)](https://join.slack.com/t/qcarchive/shared_invite/enQtNDIzNTQ2OTExODk0LTE3MWI0YzBjNzVhNzczNDM0ZTA5MmQ1ODcxYTc0YTA1ZDQ2MTk1NDhlMjhjMmQ0YWYwOGMzYzJkZTM2NDlmOGM) |
| **Foundation** | [![license](https://img.shields.io/github/license/qcdb/qcdb.svg)](https://opensource.org/licenses/BSD-3-Clause) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS%2C%20Windows%20WSL-orange.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) [![python](https://img.shields.io/badge/python-3.6+-blue.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) |

# Example

A simple example of QCDB's capabilities is as follows:

```python3
>>> import qcdb

>>> mol = qcdb.set_molecule("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

>>> qcdb.set_keywords({
    "freeze_core": True,
})
```

These input specifications plus a model chemistry can be executed with `energy()`, `gradient()`, or `hessian()` API functions along with a program specifier:

```python3
>>> ene = qcdb.energy("c4-mp2/cc-pvdz")
```
```python3
>>> ene
-76.22847548367803
>>> print(qcdb.print_variables())

  Variable Map:
  ----------------------------------------------------------------------------
  "CURRENT CORRELATION ENERGY"           =>      -0.201821821791 [Eh]
  "CURRENT ENERGY"                       =>     -76.228475483678 [Eh]
  "CURRENT REFERENCE ENERGY"             =>     -76.026653661887 [Eh]
  "HF TOTAL ENERGY"                      =>     -76.026653661887 [Eh]
  "MP2 CORRELATION ENERGY"               =>      -0.201821821791 [Eh]
  "MP2 DOUBLES ENERGY"                   =>      -0.201821821791 [Eh]
  "MP2 OPPOSITE-SPIN CORRELATION ENERGY" =>      -0.151079672317 [Eh]
  "MP2 SAME-SPIN CORRELATION ENERGY"     =>      -0.050742149474 [Eh]
  "MP2 SINGLES ENERGY"                   =>       0.000000000000 [Eh]
  "MP2 TOTAL ENERGY"                     =>     -76.228475483678 [Eh]
  "N ALPHA ELECTRONS"                    =>       5.000000000000 []
  "N ATOMS"                              =>       3.000000000000 [Eh]
  "N BASIS FUNCTIONS"                    =>      24.000000000000 []
  "N BETA ELECTRONS"                     =>       5.000000000000 []
  "N MOLECULAR ORBITALS"                 =>      24.000000000000 []
  "NUCLEAR REPULSION ENERGY"             =>       9.168193296400 [Eh]
  "SCS(N)-MP2 CORRELATION ENERGY"        =>      -0.089306183074 [Eh]
  "SCS(N)-MP2 TOTAL ENERGY"              =>     -76.115959844961 [Eh]
  "SCS-MP2 CORRELATION ENERGY"           =>      -0.198209656605 [Eh]
  "SCS-MP2 TOTAL ENERGY"                 =>     -76.224863318492 [Eh]
  ```

For a fuller record of the computation, request the [`AtomicResult`](https://molssi.github.io/QCEngine/single_compute.html?highlight=atomicresult#qcelemental.models.AtomicResult)-like return

```python3
>>> ene, ret = qcdb.energy("c4-mp2/cc-pvdz", return_wfn=True)
```
```python3
>>> ret["return_result"]
-76.228475483678

>>> ret["properties"]["mp2_correlation_energy"]
-0.201821821791

>>> ret["provenance"]["cpu"]
'Intel(R) Xeon(R) CPU E5-2630 v4 @ 2.20GHz'
```

The accustomed output (and generated input) are also accessible

```python3
>>> >>> pprint.pprint(ret["stdout"], width=200)
(' --invoking executable--\n'
 '/home/psilocaluser/gits/cfour-public/bin/xjoda\n'
 '\n'
 '\n'
 '   *************************************************************************\n'
 '         <<<     CCCCCC     CCCCCC   |||     CCCCCC     CCCCCC   >>>\n'
 '       <<<      CCC        CCC       |||    CCC        CCC         >>>\n'
 '      <<<      CCC        CCC        |||   CCC        CCC            >>>\n'
 '    <<<        CCC        CCC        |||   CCC        CCC              >>>\n'
 '      <<<      CCC        CCC        |||   CCC        CCC            >>>\n'
 '       <<<      CCC        CCC       |||    CCC        CCC         >>>\n'
 '         <<<     CCCCCC     CCCCCC   |||     CCCCCC     CCCCCC   >>>\n'
 '   *************************************************************************\n'
 '\n'
 '     ****************************************************************\n'
 '     * CFOUR Coupled-Cluster techniques for Computational Chemistry *\n'
 '     ****************************************************************\n'
 ' \n'
 '\n'
 '   Department of Chemistry                Institut fuer Physikalische Chemie\n'
 '   University of Florida                  Universitaet Mainz\n'
 '   Gainesville, FL 32611, USA             D-55099 Mainz, Germany\n'
 '\n'
 '   Department of Chemistry                Fakultaet fuer Chemie und Biowiss.\n'
 '   Johns Hopkins University               Karlsruher Institut fuer Technologie\n'
 '   Baltimore, MD 21218, USA               D-76131 Karlsruhe, Germany\n'
 '\n'
 '   Department of Chemistry                Department of Physical Chemistry\n'
 '   Southern Methodist University          Eotvos Lorand University\n'
 '   Dallas, TX 75275, USA                  H-1053 Budapest, Hungary\n'
 '\n'
 ' \n'
 '                       Version 2.1\n'
 ' \n'
 '                     psinet                                            \n'
 '                     Thu Sep  9 15:52:23 EDT 2021                      \n'
 '                     integer*8 version is running\n'
 ' \n'
 '********************************************************************************\n'
 '*                          Input from ZMAT file                                *\n'
 '********************************************************************************\n'
 'auto-generated by QCElemental from molecule H2O                                 \n'
 'O                     0.000000000000     0.000000000000    -0.124297814080      \n'
 'H                     0.000000000000    -1.434419274846     0.986348254917      \n'
 'H                     0.000000000000     1.434419274846     0.986348254917      \n'
 '                                                                                \n'
 '                                                                                \n'
 '*CFOUR(BASIS=SPECIAL                                                            \n'
 'CALC_LEVEL=MP2                                                                  \n'
 'CHARGE=0                                                                        \n'
 'COORDINATES=CARTESIAN                                                           \n'
 'DERIV_LEVEL=ZERO                                                                \n'
 'FROZEN_CORE=1                                                                   \n'
 'MEMORY_SIZE=2694957760                                                          \n'
 'MEM_UNIT=INTEGERWORDS                                                           \n'
 'MULTIPLICITY=1                                                                  \n'
 'SCF_DAMPING=0                                                                   \n'
 'SCF_MAXCYC=100                                                                  \n'
 'SPHERICAL=1                                                                     \n'
 'UNITS=BOHR)                                                                     \n'
 '                                                                                \n'
 'O:CD_1                                                                          \n'
 'H:CD_2                                                                          \n'
 'H:CD_3                                                                          \n'
 '                                                                                \n'
 '********************************************************************************\n'
 
 ...
 
 '********************************************************************************\n'
 '   The full molecular point group is C2v .\n'
 '   The largest Abelian subgroup of the full molecular point group is C2v .\n'
 '   The computational point group is C2v .\n'
 '********************************************************************************\n'
 '\n'
 '\n'
 ' ----------------------------------------------------------------\n'
 '         Coordinates used in calculation (QCOMP) \n'
 ' ----------------------------------------------------------------\n'
 ' Z-matrix   Atomic            Coordinates (in bohr)\n'
 '  Symbol    Number           X              Y              Z\n'
 ' ----------------------------------------------------------------\n'
 '     O         8         0.00000000     0.00000000    -0.12429781\n'
 '     H         1         0.00000000    -1.43441927     0.98634825\n'
 '     H         1         0.00000000     1.43441927     0.98634825\n'
 ' ----------------------------------------------------------------\n'
 ' \n'
 '   Interatomic distance matrix (Angstroms) \n'
 ' \n'
 '                 O             H             H    \n'
 '                [ 1]        [ 2]        [ 3]\n'
 '  O    [ 1]     0.00000\n'
 '  H    [ 2]     0.96000     0.00000\n'
 '  H    [ 3]     0.96000     1.51812     0.00000\n'
 ' rotcon2\n'
 '   Rotational constants (in cm-1): \n'
 '          9.4721706374            27.2629827680            14.5153365101\n'
 '   Rotational constants (in MHz): \n'
 '     283968.5715813829        817323.7761473177        435158.9020696688\n'
 '  ECPDATA file not present.   Using default ECPDATA. \n'
 '  There is   1 frozen-core orbital.\n'
 '  There are   24 basis functions.\n'
 
 ...
 
 ' --invoking executable--\n'
 '/home/psilocaluser/gits/cfour-public/bin/xvcc\n'
 ' @GETMEM-I,  Allocated   20560 MB of main memory.\n'
 '   MBPT(2) energy will be calculated.\n'
 '   The total correlation energy is -0.201821821791 a.u.\n'
 '  -----------------------------------------------------------\n'
 '       Correction         Increment           Cumulative\n'
 '  -----------------------------------------------------------\n'
 '        D-MBPT(2)       -0.201821821791     -76.228475483678\n'
 '  -----------------------------------------------------------\n'
 '  Total MBPT(2)       energy:     -76.228475483678 a.u.\n'
 ' @CHECKOUT-I, Total execution time (CPU/WALL):        0.00/       0.00 seconds.\n'
 '--executable xvcc finished with status     0 in        0.02 seconds (walltime).\n'
 '  The final electronic energy is       -76.228475483677684 a.u. \n'
 '  This computation required                            0.33 seconds (walltime).\n')
 ```

Some tutorial tests
* [test_tu1_ene.py](qcdb/tests/test_tu1_ene.py)
* [test_tu2_uhf.py](qcdb/tests/test_tu2_uhf.py)
* [test_tu3_opt.py](qcdb/tests/test_tu3_opt.py)
* [test_tu4_freq.py](qcdb/tests/test_tu4_freq.py)
* [test_tu6_cp.py](qcdb/tests/test_tu6_cp.py)

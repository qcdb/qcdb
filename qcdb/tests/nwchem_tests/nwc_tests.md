**QA tests for NWChem capabilities**

_Last update: 04/23/2020_

All tests have the prefix 'test\_' so all titles here describe the method and run.

-----   Passed  -----

Energy
- [x] dft\_xc
- [x] n2\_ccsd\_t\_
- [x] opt\_rhf\_scf
- [x] rhf\_ccsdtq (It's best not to run this one when running full test suite of nwc due to runtime)
- [x] rhf\_ccsd\_pr\_t\_
- [x] rohf\_scf\_nh2
- [x] sp\_dft\_ccsdt
- [x] sp\_rhf\_ccs\_t-2
- [x] sp-rhf-mp2-energy
- [x] sp-rhf-scf
- [x] sp\_rhf\_ccsd\_t-2
- [x] mcscf\_ch2
- [x] tce-ci-h2o
- [x] tce\_ccd
- [x] tce\_ccsd\_pr\_br\_t
- [x] tce\_mpn\_h2o
- [x] tce\_qcisd
- [x] tce\_rhf\_ccsd
- [x] tce\_uhf\_uccsd
- [x] uhf\_scf
- [x] xe\_zora

Gradient
- [x] sp\_mp2\_grad

-----   Incomplete  -----

Energy
- [ ] sp-uhf-mp2: SCS not working
- [ ] tddft_n2+

Gradient
- [ ] rhf\_scf\_grad: Gradient alignment issue
 
Properties
- [ ] tce\_ccsd\_dipole: Harvesting fix

-----   Testing  -----

Energy
- [ ] dft\_xc\_multi

Properties
- [ ] dft\_dipole

-----   Needed   -----
- [ ] direct\_mp2
- [ ] nmr\_calculations
- [ ] hessian

----- Extra information -----

-----   Update Log  -----

- 04/23/20: Create testing file, put in nwc testing as established

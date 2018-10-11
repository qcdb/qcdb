# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import absolute_import
from __future__ import print_function
import re
import sys
import struct
from collections import defaultdict
from decimal import Decimal

from ..pdict import PreservingDict 
from ..periodictable import *
from ..physconst import *
from ..exceptions import *
from ..molecule import Molecule
#from .orient import OrientMols
#from .options import conv_float2negexp
#from . import nwc_movecs

def harvest_output(outtext):
    """Function to separate portions of a NWChem output file *outtext*,
    divided by " Line search:".

    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []

    for outpass in re.split(r' Line search:',outtext,re.MULTILINE):
        psivar, nwcoord, nwgrad, version, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(nwcoord)
        pass_grad.append(nwgrad)
    
#       print ('\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n')
#       print (outpass)
#       print (psivar, nwcoord, nwgrad)
#       print (psivar, nwgrad)
#       print ('\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n')

    retindx = -1 if pass_coord[-1] else -2 

#    print ('    <<<  NW PSIVAR  >>>')
#    for item in pass_psivar[retindx]:
#       print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
#    print ('    <<<  NW COORD   >>>')
#    for item in pass_coord[retindx]:
#    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    print ('    <<<   NW GRAD   >>>')
#    for item in pass_grad[retindx]:
#       print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx], version, error


def harvest_outfile_pass(outtext):
    """Function to read NWChem output file *outtext* and parse important 
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ''
    error = ''

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Process version
    mobj = re.search(
        r'^\s+' + r'Northwest Computational Chemistry Package (NWChem)' + r'\s+' + r'(?:<version>\d+.\d+)' + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched version')
        version = mobj.group('version')

    #Process SCF 
    #1)Fail to converge
    mobj = re.search(
        r'^\s+' + r'(?:Calculation failed to converge)'+r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('failed to converge')

    #2)Calculation converged
    else:
        mobj = re.search(
            r'^\s+' + r'(?:Total SCF energy)' + r'\s+=\s*' + NUMBER +r's*$'
            ,outtext, re.MULTILINE) 
        if mobj:
            print('matched HF') 
            psivar['HF TOTAL ENERGY'] = mobj.group(1) 

    #Process Effective nuclear repulsion energy (a.u.) 
        mobj = re.search( 
            r'^\s+' + r'Effective nuclear repulsion energy \(a\.u\.\)' + r'\s+' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched NRE')
            #print (mobj.group(1)) 
            psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

        #Process DFT (RDFT, RODFT,UDFT, SODFT)
        mobj = re.search(
            r'^\s+' + r'Total DFT energy' + r'\s' + NUMBER + r'\s*$',
            outtext, re.MULTILINE)
        if mobj:
            print ('matched DFT')
            #print (mobj.group(1))
            psivar ['DFT TOTAL ENERGY'] = mobj.group(1)

    #MCSCF 
        mobj = re.search(
            r'^\s+' + r'Total SCF energy' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'One-electron energy' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'Two-electron energy' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'Total MCSCF energy' + r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE) #MCSCF
        if mobj:
            print('Total MCSCF energy') #MCSCF energy calculation 
            psivar ['SCF TOTAL ENERGY'] = mobj.group(1)
            psivar ['MCSCF TOTAL ENERGY'] = mobj.group(4)

        #Process MP2 (Restricted, Unrestricted(RO n/a))
        #1)SCF-MP2
        mobj = re.search(
        r'^\s+' + r'SCF energy' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Correlation energy' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Singlet pairs' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Triplet pairs' + r'\s+' + NUMBER + r'\s*'+ 
        r'^\s+' + r'Total MP2 energy' + r'\s+' + NUMBER + r'\s*$'
        ,outtext, re.MULTILINE) #MP2
        if mobj:
            print ('matched scf-mp2')
            psivar ['CURRENT REFERENCE ENERGY'] = mobj.group(1) 
            psivar ['MP2 CORRELATION ENERGY'] = mobj.group(2)
            psivar ['MP2 TOTAL ENERGY'] = mobj.group(5) 

        mobj = re.search(
        r'^\s+' + r'Same spin pairs'+ r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Same spin scaling factor'+ r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Opposite spin pairs'+ r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Opposite spin scaling fact.'+ r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'SCS-MP2 correlation energy' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Total SCS-MP2 energy'+r'\s+' + NUMBER + r'\s*$'
        ,outtext, re.MULTILINE)  # SCS-MP2
        if mobj:
            print('matched scs-mp2')
            psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * Decimal(mobj.group(2))
            psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) * Decimal(mobj.group(4))
            psivar['SCS-MP2 CORRELATION ENERGY'] = mobj.group(5)
            psivar['SCS-MP2 TOTAL ENERGY'] = mobj.group(6)
            
            print (mobj.group(1)) #ess
            print (mobj.group(2)) #fss
            print (mobj.group(3)) #eos
            print (mobj.group(4)) #fos
            print (mobj.group(5)) #scs corl
            print (mobj.group(6)) #scs-mp2

        #2) DFT-MP2
        mobj = re.search(
            r'^\s+' + r'DFT energy' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'Unscaled MP2 energy' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'Total DFT+MP2 energy' + r'\s+' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE) 
        if mobj:
            print ('matched dft-mp2')
            psivar ['CURRENT REFERENCE ENERGY'] = mobj.group(1)
            psivar ['MP2 CORRELATION ENERGY'] = mobj.group(2)
            psivar ['MP2 TOTAL ENERGY'] = mobj.group(3)

        #3) MP2 with CCSD or CCSD(T) calculation (through CCSD(T) directive)
        mobj = re.search(
            r'^\s+' + r'MP2 Energy \(coupled cluster initial guess\)' + r'\s*'+
            r'^\s+' + r'------------------------------------------' + r'\s*'+
            r'^\s+' + r'Reference energy:' + r'\s+' + NUMBER + r'\s*'+ 
            r'^\s+' + r'MP2 Corr\. energy:' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'Total MP2 energy:' + r'\s+' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE)

        if mobj:
            print ('matched coupled cluster-mp2')
            psivar ['MP2 CORRELATION ENERGY'] = mobj.group(2) 
            psivar ['MP2 TOTAL ENERGY'] = mobj.group(3) 
 
        #Process CC calculation through tce [dertype] command
        cc_name = ''
        mobj = re.search(
            r'^\s+' + r'Iterations converged' + r'\s*' +
            r'^\s+' + r'(.*?)' + r' correlation energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*'+
            r'^\s+' + r'(.*?)' + r' total energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE) 

        if mobj:
            cc_name = mobj.group(1)
            print('matched %s' % (mobj.group(1))) 
            print (mobj.group(1)) #cc_name
            print (mobj.group(2)) #cc_name corl. energy
            print (mobj.group(4)) #cc_name total energy

            psivar['%s CORRELATION ENERGY' % (mobj.group(1))] = mobj.group(2)
            psivar['%s TOTAL ENERGY' % (mobj.group(1))] = mobj.group(4)

         # Process CC '()' correction part through tce [dertype] command
        mobj = re.search(
            r'^\s+' + cc_name + r'\('+r'(.*?)'+r'\)'+ r'\s+' + r'correction energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*'+
            r'^\s+' + cc_name + r'\('+r'(.*?)'+r'\)'+ r'\s+' + r'correlation energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*'+
            r'^\s+' + cc_name + r'\('+r'(.*?)'+r'\)'+ r'\s+' + r'total energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE)

        if mobj:
            print('matched %s(%s)' % (cc_name, mobj.group(1)))
                #print (mobj.group(1)) #correction ()
                #print (mobj.group(2)) #correction () correction energy
                #print (mobj.group(4)) #correction () corl. energy 
                #print (mobj.group(6)) #correction () total energy

            psivar['(%s) CORRECTION ENERGY'  % (mobj.group(1))] = mobj.group(2)
            psivar['%s(%s) CORRELATION ENERGY' % (cc_name, mobj.group(1))] = mobj.group(4)
            psivar['%s(%s) TOTAL ENERGY'  % (cc_name, mobj.group(1))] = mobj.group(6) 

        #Process CCSD/CCSD(T) using nwchem CCSD/CCSD(T) [dertype] command
        mobj = re.search(
        r'^\s+' + r'-----------' + r'\s*' +
        r'^\s+' + r'CCSD Energy' + r'\s*' +
        r'^\s+' + r'-----------' + r'\s*' +
        r'^\s+' + r'Reference energy:' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'CCSD corr\. energy:' + r'\s+' + NUMBER + r'\s*'+
        r'^\s+' + r'Total CCSD energy:' + r'\s+' + NUMBER + r'\s*$'
        ,outtext, re.MULTILINE | re.DOTALL)

        if mobj:
            print ('matched ccsd')
            psivar['CCSD CORRELATION ENERGY'] = mobj.group(2)
            psivar['CCSD TOTAL ENERGY'] = mobj.group(3)

        mobj = re.search(
        r'^\s+' + r'--------------' + r'\s*' +
        r'^\s+' + r'CCSD\(T\) Energy' + r'\s*' +
        r'^\s+' + r'--------------' + r'\s*' +
        r'(?:.*?)' + 
        r'^\s+' + r'\(T\) corr\. energy:' + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'Total CCSD\(T\) energy:' + r'\s+' + NUMBER + r'\s*$' 
        ,outtext, re.MULTILINE | re.DOTALL)

        if mobj:
            print ('matched ccsd(t)')
            psivar['(T) CORRECTION ENERGY'] = mobj.group(1)
            psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(2)) - psivar['HF TOTAL ENERGY'] 
            psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(2)

        #Process EOM-[cc_name] #nwchem_tce_dipole = false
        # Parsed information: each symmetry, root excitation energy in eV and total energy in hartree
        # psivar name might need to be fixed
        # each root excitation energy is extracted from the last iteration of right hand side 
        mobj = re.findall(
        r'^\s+' + r'Excited-state calculation'+r'\s+'+r'\(\s+'+ r'(.*?)'+ r'\s+'+r'symmetry\)'+ r'\s*'+
        r'^\s+' + r'========================================='+ r'\s*'+
        r'(?:.*?)' +  #Skip to line before right-hand side iterations
        r'^\s+' + r'EOM-'+ cc_name +r' right-hand side iterations' + r'\s*' +
        r'(?:.*?)'+
        r'^\s+' + r'Iteration' + r'\s+\d+\s+' + r'using' + r'\s+\d+\s+' + r'trial vectors' + r'\s*' + 
        r'((?:\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+(?:(\s+\d+\.\d+\s+\d+\.\d+)?)\s*\n)+)'+
        r'^\s+' + r'--------------------------------------------------------------' + r'\s*'+
        r'^\s+' + r'Iterations converged' + r'\s*$'
        ,outtext, re.MULTILINE | re.DOTALL)

        if mobj:
            ext_energy = {} #dic
            
            ext_energy_list = []
            for mobj_list in mobj:
                print ('matched EOM-%s - %s symmetry' %(cc_name, mobj_list[0])) #cc_name, symmetry
                print (mobj_list)
                count = 0
                for line in mobj_list[1].splitlines():
                    lline = line.split()
                    print (lline[1]) #in hartree
                    print (lline[2]) #in eV
                    count +=1

                    print ('matched excitation energy #%d - %s symmetry' %(count, mobj_list[0]))

                    ext_energy_list.append(lline[1]) #Collect all excitation energies 

                    sym = str(mobj_list[0]) 
                    ext_energy.setdefault(sym,[])
                    ext_energy[sym].append(lline[1]) #Dictionary: symmetries(key), energies(value)
  
            ext_energy_list.sort(key=float)

            for nroot in range(len(ext_energy_list)):
                for k, e_val in ext_energy.items():
                    if ext_energy_list[nroot] in e_val:
                        symm = k
                        psivar ['EOM-%s ROOT 0 -> ROOT %d EXCITATION ENERGY - %s SYMMETRY' %(cc_name, nroot+1, symm)] = \
                            ext_energy_list[nroot] #in hartree
                        psivar ['EOM-%s ROOT 0 -> ROOT %d TOTAL ENERGY - %s SYMMETRY' %(cc_name, nroot+1, symm)] = \
                            psivar['%s TOTAL ENERGY' %(cc_name)] + Decimal(ext_energy_list[nroot]) #in hartree
        #TCE_CR_EOMCCSD(T) information
        mobj = re.search(
            r'^\s+' + r'CR-' + r'(w+)' + r'\s+' + r'total energy / hartree' + r'\s+' + NUMBER + r'\s*' +  
            r'^\s+' + r'CR-' + r'(w+)' + r'\s+' + r'excitation energy \(eV\)' +r'\s+' + NUMBER + r'\s*$',
            outtext, re.MULTILINE) 
        if mobj:
            print(mobj)
            print("matched CR =")
       
        
        #TCE- ROHF and UHF
        # mobj = re.findall(
       #     r'^\s+' + 'Total SCF energy' + '^\s+' + NUMBER +
        #    r'^\s*' + r'^\s+' + r'CCSD correlation energy / hartree' + r'\s+' + NUMBER +
        #    r'^\s*' + r'\s+' + r'CCSD total energy / hartree' + r'\s+'+ NUMBER + r's*$', outtext, re.MULTILINE)
       # if mobj:
        #    print(mobj) #print list
         #   print('CCSD Energy = ')
          #  psivar['CCSD CORRELATION ENERGY'] = mobj.group(5)
           # psivar['CCSD TOTAL ENERGY'] = mobj.group(6)


                 # No symmetry
#                psivar ['EOM-%s ROOT 0 -> ROOT %d EXCITATION ENERGY' %(cc_name, nroot+1)] = \
#                    ext_energy_list[nroot] #in hartree
#                psivar ['EOM-%s ROOT 0 -> ROOT %d TOTAL ENERGY' %(cc_name, nroot+1)] = \
#                    psivar['%s TOTAL ENERGY' %(cc_name)] + Decimal(ext_energy_list[nroot]) 

        #Process TDDFT
#       1) Spin allowed
        mobj = re.findall(
        r'^\s+' + r'----------------------------------------------------------------------------'+ r'\s*'+
        r'^\s+' + r'Root'+ r'\s+'+r'(\d+)'+r'\s+'+r'(\w+)'+r'\s+'+r'(.*?)'+r'\s+'+ NUMBER + r'\s+a\.u\.\s+'+ NUMBER + r'\s+eV\s*'+
        r'^\s+' + r'----------------------------------------------------------------------------'+ r'\s*'+
        r'^\s+' + r'Transition Moments' + r'\s+X\s+'+ NUMBER + r'\s+Y\s+'+ NUMBER +r'\s+Z\s+'+ NUMBER + r'\s*'+
        r'^\s+' + r'Transition Moments' + r'\s+XX\s+'+ NUMBER + r'\s+XY\s+'+ NUMBER + r'\s+XZ\s+'+ NUMBER + r'\s*'+
        r'^\s+' + r'Transition Moments' + r'\s+YY\s+'+ NUMBER + r'\s+YZ\s+'+ NUMBER + r'\s+ZZ\s+'+ NUMBER + r'\s*$'
        ,outtext, re.MULTILINE) 

        if mobj:
            print ('matched TDDFT with transition moments')
            for mobj_list in mobj:
                #print (mobj_list) 
                #print (mobj_list[0]) #root number
                #print (mobj_list[1]) #singlet or triplet?
                #print (mobj_list[2]) #symmetry

             #### temporary psivars ####    
                psivar ['TDDFT ROOT %s %s %s EXCITATION ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = mobj_list[3]  # in a.u.
                psivar ['TDDFT ROOT %s %s %s EXCITED STATE ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = \
                    psivar ['DFT TOTAL ENERGY'] + Decimal(mobj_list[3]) 
                psivar ['TDDFT ROOT %s DIPOLE X' %(mobj_list[0])] = mobj_list[5]
                psivar ['TDDFT ROOT %s DIPOLE Y' %(mobj_list[0])] = mobj_list[6]
                psivar ['TDDFT ROOT %s DIPOLE Z' %(mobj_list[0])] = mobj_list[7]
                psivar ['TDDFT ROOT %s QUADRUPOLE XX' %(mobj_list[0])] = mobj_list[8]
                psivar ['TDDFT ROOT %s QUADRUPOLE XY' %(mobj_list[0])] = mobj_list[9]
                psivar ['TDDFT ROOT %s QUADRUPOLE XZ' %(mobj_list[0])] = mobj_list[10]
                psivar ['TDDFT ROOT %s QUADRUPOLE YY' %(mobj_list[0])] = mobj_list[11]
                psivar ['TDDFT ROOT %s QUADRUPOLE YZ' %(mobj_list[0])] = mobj_list[12]
                psivar ['TDDFT ROOT %s QUADRUPOLE ZZ' %(mobj_list[0])] = mobj_list[13]

#       2) Spin forbidden
        mobj = re.findall(
        r'^\s+' + r'----------------------------------------------------------------------------'+ r'\s*'+
        r'^\s+' + r'Root'+ r'\s+'+r'(\d+)'+r'\s+'+r'(\w+)'+r'\s+'+r'(.*?)'+r'\s+'+ NUMBER + r'\s+a\.u\.\s+'+ NUMBER + r'\s+eV\s*'+
        r'^\s+' + r'----------------------------------------------------------------------------'+ r'\s*'+
        r'^\s+' + r'Transition Moments' + r'\s+Spin forbidden\s*$'
        ,outtext, re.MULTILINE)

        if mobj:
            print ('matched TDDFT - spin forbidden')
            for mobj_list in mobj:
             #### temporary psivars ####    
                psivar ['TDDFT ROOT %s %s %s EXCITATION ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = mobj_list[3]  # in a.u.
                psivar ['TDDFT ROOT %s %s %s EXCITED STATE ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = \
                    psivar ['DFT TOTAL ENERGY'] + Decimal(mobj_list[3])
            if mobj:
                print('Non-variation initial energy') #prints out energy, 5 counts
    
        #Process geometry
        # 1) CHARGE
        # Read charge from SCF module
        mobj = re.search(
            r'^\s+' + r'charge          =' + r'\s+'+ NUMBER + r'\s*$'
            ,outtext,re.MULTILINE | re.IGNORECASE)

        if mobj:
            print ('matched charge')
            out_charge = int(float(mobj.group(1)))

        # Read charge from General information (not scf module)
        mobj = re.search(
        r'^\s+' + r'Charge           :' + r'\s+'+ r'(-?\d+)'+ r'\s*$' 
        ,outtext,re.MULTILINE | re.IGNORECASE)

        if mobj:
            print ('matched charge')
            out_charge = int(float(mobj.group(1)))
        
        # 2) MULTIPLICITY
        # Read multiplicity from SCF module
        mobj = re.search(
        r'^\s+' + r'open shells     =' + r'\s+'+ r'(\d+)'+ r'\s*$' 
        ,outtext,re.MULTILINE | re.IGNORECASE)

        if mobj: 
            print ('matched multiplicity')
            out_mult = int(mobj.group(1)) + 1

        # Read multiplicity from SCF module through alpha, beta electrons
        mobj = re.search(
        r'^\s+' + r'alpha electrons =' + r'\s+'+ r'(\d+)'+ r'\s*'+
        r'^\s+' + r'beta  electrons =' + r'\s+'+ r'(\d+)'+ r'\s*$'
        ,outtext,re.MULTILINE | re.IGNORECASE)   

        if mobj:
            print ('matched multiplicity')
            out_mult = int(mobj.group(1))-int(mobj.group(2)) + 1 #nopen + 1
  
        # Read multiplicity from General information (not scf module)
        mobj = re.search(
        r'^\s+' + r'Spin multiplicity:' + r'\s+'+ r'(\d+)' + r'\s*$'
        ,outtext,re.MULTILINE | re.IGNORECASE)

        if mobj:
            print ('matched multiplicity')
            out_mult = int(mobj.group(1)) 

        #3) Initial geometry   
        mobj = re.search(
        r'^\s+' + r'Geometry' + r'.*' + r'\s*' +
        r'^\s+' + r'(?:-+)\s*'+
        r'\s+' + r'\n' +
        r'^\s' + r'Output coordinates in ' + r'(.*?)' + r'\s' + r'\(scale by'+ r'.*' + r'\s'+ r'to convert to a\.u\.\)'+
        r'\s+' + r'\n' +
        r'^\s+' + r'No\.\       Tag          Charge          X              Y              Z'+ r'\s*'+
        r'^\s+' + r'---- ---------------- ---------- -------------- -------------- --------------'+r'\s*'+
        r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-z]*)+\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'+r'\s*$'
        ,outtext,re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched geom')
            print(mobj.groups())

            #dinky molecule w/ charge and multiplicity
            if mobj.group(1) == 'angstroms':
                #print ('unit is angstroms')
                molxyz = '%d \n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult) # unit = angstrom
                for line in mobj.group(2).splitlines():
                    lline = line.split()
                    molxyz += '%s %16s %16s %16s\n' %(lline[-5], lline[-3], lline[-2], lline[-1])
                    # Jiyoung was collecting charge (-4)? see if this is ok for ghosts
                                                  # Tag    ,    X,        Y,        Z
                #psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)
                psivar_coord = Molecule.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)

        
            else: # unit = a.u.
                #print ('unit is au') 
                molxyz = '%d au\n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult)
                for line in mobj.group(2).splitlines():
                    lline = line.split()
                    molxyz += '%s %16s %16s %16s\n' %(lline[-4], lline[-3], lline[-2], lline[-1])
                                                      # Tag    ,    X,        Y,        Z
                #psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)
                psivar_coord = Molecule.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)
            
        #Process gradient
        mobj = re.search(
            r'^\s+' + r'.*'+ r'ENERGY GRADIENTS' + r'\s*'+
            r'\s+' + r'\n'+
            r'^\s+' + r'atom               coordinates                        gradient' + r'\s*'+
            r'^\s+' + r'x          y          z           x          y          z' + r'\s*' +
            r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-x]*)+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'+r'\s*$'
            ,outtext, re.MULTILINE)

        if mobj:
            print ('matched molgrad')
            atoms = []
            psivar_grad = []
            for line in mobj.group(1).splitlines():
                lline = line.split()  # Num, Tag, coord x, coord y, coord z, grad x, grad y, grad z
                #print (lline)
                if lline == []:
                    pass
                else: 
                    atoms.append(lline[1]) #Tag
                    psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

        #Process dipole 
        mobj = re.search(
            r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s' + r'A\.U\.'+ r'\s*' + 
            r'^\s+' + r'DMX' + r'\s+' + NUMBER + r'.*' + r'\s*' +
            r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' + r'\s*' +
            r'^\s+' + r'DMZ' + r'\s+' + NUMBER + r'.' + r'\s*' + 
            r'^\s+' + r'.*' + r'\s*' + 
            r'^\s+' + r'Total dipole' + r'\s+' + NUMBER + r'\s' + r'A\.U\.'+ r'\s*'+
            r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s' + r'Debye\(s\)' + r'\s*'+
            r'^\s+' + r'DMX' + r'\s+' + NUMBER + r'.*' + r'\s*' +
            r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' + r'\s*' +
                r'^\s+' + r'DMZ' + r'\s+' + NUMBER + r'.*' + r'\s*' +
            r'^\s+' + r'.*' + r'\s*' +
            r'^\s+' + r'Total dipole' + r'\s+' + NUMBER + r'\s' + r'DEBYE\(S\)' + r'\s*$'
            ,outtext,re.MULTILINE)
                
        if mobj:
            print ('matched total dipole')  
           #print (mobj.group(5), 'a.u.') 
           #print (mobj.group(10),'debye(s)')

           #UNIT = DEBYE(S)
            psivar['CURRENT DIPOLE X'] = mobj.group(7)
            psivar['CURRENT DIPOLE Y'] = mobj.group(8)
            psivar['CURRENT DIPOLE Z'] = mobj.group(9)
           # total?             

        #Process error code
            mobj = re.search(
                r'^\s+' + r'current input line \:' +r'\s*'+
                r'^\s+'+ r'([1-9][0-9]*)' + r'\:' +r'\s+' + r'(.*)'+ r'\s*'+
                r'^\s+'r'------------------------------------------------------------------------' + r'\s*'+
                r'^\s+'r'------------------------------------------------------------------------' + r'\s*'+
                r'^\s+' + r'There is an error in the input file' + r'\s*$'
                ,outtext,re.MULTILINE)
            if mobj:
                print('matched error')
               #print (mobj.group(1)) #error line number
           #print (mobj.group(2)) #error reason 
            psivar['NWCHEM ERROR CODE'] = mobj.group(1)    
            # TODO process errors into error var

    # Process CURRENT energies (TODO: needs better way)
    if 'HF TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['HF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['HF TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar and 'MP2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']

    if 'DFT TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['DFT TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['DFT TOTAL ENERGY']

    # Process TCE CURRENT energies
      # Need to be fixed  
      # HOW TO KNOW options['NWCHEM']['NWCHEM_TCE']['value']?
      # TODO: CURRENT ENERGY = TCE ENERGY
    if ('%s TOTAL ENERGY' % (cc_name) in psivar and \
       ('%s CORRELATION ENERGY' % (cc_name) in psivar)):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (cc_name)]
        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (cc_name)]
  
    if 'CCSD(T) TOTAL ENERGY' in psivar and 'CCSD(T) CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']
    
    if 'CISD CORRELATION ENERGY' in psivar:
        psivar['CISD CORRELATION ENERGY'] = psivar['CISD CORRELATION ENERGY']
    
    if 'CISD TOTAL ENERGY' in psivar:
        psivar['CISD TOTAL ENERGY'] = psivar['CISD TOTAL ENERGY']
   
    if 'CISDT CORRELATION ENERGY' in psivar:
        psivar['CISDT CORRELATION ENERGY'] = psivar['CISDT CORRELATION ENERGY']
               
    if 'CISDT TOTAL ENERGY' in psivar:
        psivar['CISDT TOTAL ENERGY'] = psivar['CISDT TOTAL ENERGY']

    if 'MP2 CORRELATION ENERGY' in psivar:
        psivar['MP2 CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar:
        psivar['MP2 TOTAL ENERGY'] = psivar['MP2 TOTAL ENERGY']

    if 'MP3 CORRELATION ENERGY' in psivar:
        psivar['MP3 CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY']

    if 'MP3 TOTAL ENERGY' in psivar:
        psivar['MP3 TOTAL ENERGY'] = psivar['MP3 TOTAL ENERGY']

    if 'MP4 CORRELATION ENERGY' in psivar:
        psivar['MP4 CORRELATION ENERGY'] = psivar['MP4 CORRELATION ENERGY']

    if 'MP4 TOTAL ENERGY' in psivar:
        psivar['MP4 TOTAL ENERGY'] = psivar['MP4 TOTAL ENERGY']

    if ('EOM-%s TOTAL ENERGY' % (cc_name) in psivar) and \
       ('%s EXCITATION ENERGY' %(cc_name) in psivar):
        psivar['CURRENT ENERGY'] = psivar['EOM-%s TOTAL ENERGY' %(cc_name)]
        psivar['CURRENT EXCITATION ENERGY'] = psivar['%s EXCITATION ENERGY' %(cc_name)] 

    return psivar, psivar_coord, psivar_grad, version, error

def harvest_hessian(hess):
    """Parses the contents *hessian* of the NWChem hess file into a hessian array.
    Hess file name has to be "nwchem.hess". (default)

    """
    hess = hess.splitlines() 
        

def muster_inherited_options(ropts, verbose=1):
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> cfour/memory_size [MB]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = [int(0.000001 * qopt.value), 'mb']
        print('\n\nMEMORY', mem, '\n\n')
        ropts.suggest('NWCHEM', 'MEMORY', mem, **kwgs)


def muster_memory(mem):
    """Transform input *mem* in MB into nwchem-speak 
    
    """
    text = ''

    options = defaultdict(lambda: defaultdict(dict))
    options['NWCHEM']['NWCHEM_MEMORY']['value'] = [int(mem), 'mb'] #Total memory set in mb

    for item in options['NWCHEM']:
        options['NWCHEM'][item]['clobber'] = True

    return text, options
    
#May need to adjust options for qcdb driver ATL
def muster_psi4options(opt, mol):
    """Translate psi4 keywords *opt* that have been explicitly set into
    their NWChem counterparts. Since explicitly set NWChem module keyword
    values will always be used preferentially to these inferred from
    psi4, the 'clobber' property is set to False.
    *mol* = qcdb.Molecule(molecule.create_psi4_string_from_molecule())   

    """
    text = ''
    options = defaultdict(lambda:defaultdict(dict))

    mult = mol.multiplicity()
   #print ('mult: ', mult)  
        
    if 'SCF' in opt:
        if 'REFERENCE' in opt['SCF']:
            if opt['SCF']['REFERENCE']['value'] == 'UKS':
                options['NWCHEM']['NWCHEM_DFT']['value'] = 'ODFT'

            elif opt['SCF']['REFERENCE']['value'] == 'RKS': 
                if mult != 1:   
                    options['NWCHEM']['NWCHEM_DFT']['value'] = 'RODFT'
                else:
                    pass
            else:
                options['NWCHEM']['NWCHEM_SCF']['value'] = \
                        {'RHF': 'RHF',
                         'UHF': 'UHF',
                         'ROHF': 'ROHF'}[opt['SCF']['REFERENCE']['value']]

        if 'SCF_TYPE' in opt['SCF']:
            if opt['SCF']['SCF_TYPE']['value'] == 'DIRECT':
                options['NWCHEM']['NWCHEM_DFT_DIRECT']['value'] = 'TRUE'
                options['NWCHEM']['NWCHEM_SCF_DIRECT']['value'] = 'TRUE'
         
        if 'D_CONVERGENCE' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] = \
                ['density', opt['SCF']['D_CONVERGENCE']['value']]
            #This is not correct so has to be fixed
            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
                opt['SCF']['D_CONVERGENCE']['value']

        if 'E_CONVERGENCE' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] += \
                ['energy', opt['SCF']['E_CONVERGENCE']['value']]
            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
                opt['SCF']['E_CONVERGENCE']['value']

        if 'MAXITER' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_ITERATIONS']['value'] = \
               opt['SCF']['MAXITER']['value']   
            options['NWCHEM']['NWCHEM_SCF_MAXITER']['value'] = \
               opt['SCF']['MAXITER']['value']

        if 'DFT_FUNCTIONAL' in opt['SCF']:
            dft_dispersion_list = \
                ['B3LYP-D3','BLYP-D3','BP86-D3','M05-D3','M05-2X-D3', \
                'PBE0-D3']
            if opt['SCF']['DFT_FUNCTIONAL']['value'] in dft_dispersion_list:
                options['NWCHEM']['NWCHEM_DFT_XC']['value'], options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = \
            muster_dft_w_dispersion(opt)

            else: 
                options['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(opt)
  
        if 'DFT_RADIAL_POINTS' and 'DFT_SPHERICAL_POINTS' in opt['SCF']:
            rad_val = opt['SCF']['DFT_RADIAL_POINTS']['value']
            sph_val = opt['SCF']['DFT_SPHERICAL_POINTS']['value']
            
            sph_val_list = \
                [38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354, \
                2702,3074,3470,3890,4334,4802,5294,5810]
            
            iangquad = sph_val_list.index(sph_val) + 1

            grid_val = []
            grid_val.append(rad_val)
            grid_val.append(iangquad)

            #print ('grid_val? ', grid_val)

            options['NWCHEM']['NWCHEM_DFT_GRID']['value'] = grid_val            
            

    for item in options['NWCHEM']:
        options['NWCHEM'][item]['clobber'] = False
    return text, options

def muster_dft_functionals(opt):
    """Translate psi4 dft functionals w/o dispersion into NWChem counterparts. 
    
    """

    text = ''
    options = defaultdict(lambda:defaultdict(dict))
    dft_functional = ''

    dft_same_name_list = ['B3LYP','PBE0','M05','M05-2X',\
      'FT97','HCTH','HCTH120', 'DLDF', \
      'HCTHP14','HCTH407P','B1B95','M06','M06-2X','M06-HF',\
      'M08-HX','M08-SO','M11','M11-L','MPW1B95','MPW1K','MPWB1K',\
      'PW6B95','PWB6K', 'BOP']
   
    val = opt['SCF']['DFT_FUNCTIONAL']['value'] 
    if val in dft_same_name_list:
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
                    {'B3LYP': 'B3LYP',
                     'PBE0' : 'PBE0',
                       'M05': 'M05',
                    'M05-2X':'M05-2X',
                     'FT97' : 'FT97',
                     'HCTH' : 'HCTH',
                   'HCTH120': 'HCTH120',
                     'DLDF' : 'DLDF',
                   'HCTHP14': 'HCTHP14',
                  'HCTH407P': 'HCTH407P',
                     'B1B95': 'B1B95',
                       'M06': 'M06',
                    'M06-2X': 'M06-2X',
                    'M06-HF': 'M06-HF',
                    'M08-HX': 'M08-HX',
                    'M08-SO': 'M08-SO',
                       'M11': 'M11',
                     'M11-L': 'M11-L',
                   'MPW1B95': 'MPW1B95',
                    'MPW1K' : 'MPW1K',
                    'MPWB1K': 'MPWB1K',
                    'PW6B95': 'PW6B95',
                     'PWB6K': 'PWB6K',
                      'BOP' : 'BOP'
                    }[opt['SCF']['DFT_FUNCTIONAL']['value']]

    elif (val =='HF'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 1.0]

    elif (val == 'B2PLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.53, 'becke88', 0.47, 'lyp', 0.73, 'mp2', 0.27]

    elif (val == 'BLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'lyp']
    
    elif (val == 'BP86'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'perdew86']

    elif (val == 'PBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xpbe96', 'cpbe96']

    elif (val=='PW91'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xperdew91', 'perdew91']

    elif (val == 'B97-1'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-1']
    elif (val == 'B97-2'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-2']
    elif (val == 'B97-GGA1'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1']

    elif (val == 'TPSSH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xctpssh']
    
    elif (val == 'BHANDH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['beckehandh']
    
    elif (val =='CAM-B3LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
        ['xcamb88','1.00','lyp','0.81','vwn_5','0.19','hfexch','1.00\n','cam','0.33',\
        'cam_alpha','0.19','cam_beta','0.46']

    elif (val =='B1LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.25,'becke88', 0.75, 'lyp']

    elif (val =='B1PW91'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.25,'becke88', 0.75, 'perdew91']

    elif (val == 'B3LYP5'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['vwn_5', 0.19, 'lyp', 0.81, 'HFexch', 0.20, 'slater',0.8,'becke88','nonlocal',0.72] 

    elif (val == 'B86BPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke86b', 'cpbe96']

    elif (val == 'B97') or (val == 'B97-0'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97','HFexch', 0.1943]        
    elif (val == 'B97-1P'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1','HFexch', 0.1500]
    elif (val == 'B97-3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-3','HFexch', 0.2693]

    elif (val == 'BHANDHLYP') or (val == 'BHHLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke88',0.500,'HFexch', 0.500,'lyp']

    elif (val == 'LRC-WPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',1.00,'cpbe96',1.0,'hfexch','1.00\n',\
          'cam',0.3,'cam_alpha',0.00,'cam_beta',1.00]
    elif (val == 'LRC-WPBEH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',0.80,'cpbe96',1.0,'hfexch','1.00\n',\
          'cam',0.2,'cam_alpha',0.20,'cam_beta',0.80]

    elif (val == 'MPW1PW'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91',0.75,'HFexch', 0.25,'perdew91']
    elif (val == 'MPWLYP1M'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91',0.95,'HFexch', 0.05,'lyp']
    elif (val == 'MPWLYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91','vwn_5', 0.12,'lyp',0.88]

    elif (val == 'PBE0-13'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',0.6667,'HFexch',0.3333,'cpbe96']
    elif (val == 'PBEH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',0.75,'HFexch',0.25,'cpbe96']
    elif (val == 'PBELYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96','vwn_5', 0.26,'lyp',0.74]

    elif (val == 'PW86PBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xperdew86','cpbe96']

    elif (val == 'TPSSLYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xtpss03','vwn_5',0.26,'lyp',0.74]

    elif (val == 'XLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
        ['slater','-0.0690','becke88',0.722,'xperdew91',0.347,'lyp']

    elif (val == 'WPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0,\
             '\ncam', 0.40, 'cam_alpha', 0.0, 'cam_beta', 1.0]

    elif (val == 'SB98-1A'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke98','HFexch',0.229]
 
    elif (val == 'PBEH3C'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',0.58,'HFexch', 0.42,'cpbe96']
    elif (val == 'PBE1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96','vwn_5', 0.26,'cpbe96',0.74]
    elif (val == 'PBE0-2'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',0.2063,'HFexch',0.7937,'cpbe96',0.5,'mp2',0.5]
    elif (val == 'O3LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['slater',0.0706,'optx',0.8133,'HFexch',0.1161,'vwn_5', 0.19,'lyp',0.81]
    elif (val == 'HSE03'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96','1.0\n',\
          'cam',0.1061,'cam_alpha',0.0,'cam_beta',0.25]
    elif (val == 'HSE06'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96',1.0,'srhfexch','0.25\n',\
          'cam',0.11,'cam_alpha',0.0,'cam_beta',1.0]    
    elif (val == 'WPBE0'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0, \
            '\ncam', 0.30, 'cam_alpha', 0.25, 'cam_beta', 0.75] 

    else:
        val = str(val)
        raise ValidationError("""Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """ % (val))

    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']

    return dft_functional

def muster_dft_w_dispersion(opt):
    """Translate psi4 dft functionals w/dispersion into NWChem counterparts. 
    
    """

    text = ''
    options = defaultdict(lambda:defaultdict(dict))
    dft_functional = ''
    dft_dispersion = ''

    val = opt['SCF']['DFT_FUNCTIONAL']['value']

    if (val == 'B3LYP-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['b3lyp']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'BLYP-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'lyp']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'BP86-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'perdew86']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'M05-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'M05-2X-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05-2x']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
  
    elif (val == 'PBE0-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['pbe0']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    else:
        val = str(val)
        raise ValidationError("""Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """ % (val))

    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']
    dft_dispersion = options['NWCHEM']['NWCHEM_DFT_DISP']['value']

    return dft_functional, dft_dispersion

def muster_modelchem(name, dertype):
    """Transform calculation method *name* and derivative level *dertype*
    into options for NWChem. While deliberately requested pieces,
    generally nwchem__nwchem_task ...  is set to complain if contradicted 
    ('clobber' set to True), 
    other 'recommended' settings, can be countermanded by keywords in input file
    ('clobber' set to False). Occasionally, we want these pieces to actually
    overcome keywords in input file ('superclobber' set to True). 

    """

    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))
    
    if lowername == 'nwchem':
        pass
    elif lowername == 'nwc-scf':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'energy'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'gradient'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'hessian'
    
    elif lowername == 'nwc-hf':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'energy'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'gradient'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_SCF']['value'] = 'hessian'
            
    elif lowername == 'nwc-dft':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'energy'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'gradient'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'hessian'

    elif lowername == 'nwc-mp2':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_MP2']['value'] = 'energy'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_MP2']['value'] = 'gradient'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_MP2']['value'] = 'hessian'

    elif lowername == 'nwc-ccsd':
        if dertype == 0:
            options['NWCHEM']['NWCHEM_TASK_CCSD']['value'] = 'energy'
        elif dertype == 1:
            options['NWCHEM']['NWCHEM_TASK_CCSD']['value'] = 'gradient'
        elif dertype == 2:
            options['NWCHEM']['NWCHEM_TASK_CCSD']['value'] = 'hessian'

    elif lowername == 'nwc-ccsd(t)':
        if dertype == 0:
            options['NWCHEM']['NWCHEM_TASK_CCSD(T)']['value'] = 'energy'
        elif dertype == 1:
            options['NWCHEM']['NWCHEM_TASK_CCSD(T)']['value'] = 'gradient'
        elif dertype == 2:
            options['NWCHEM']['NWCHEM_TASK_CCSD(T)']['value'] = 'hessian'

    elif lowername == 'nwc-ccsdt':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'energy'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'gradient'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'hessian'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'

    elif lowername == 'nwc-ccsdtq':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'energy'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'gradient'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'hessian'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'

    elif lowername == 'nwc-eom-ccsd':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'energy'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsd'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'gradient'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsd'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'hessian'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsd'

    elif lowername == 'nwc-eom-ccsdt':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'energy'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'gradient'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'hessian'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdt'

    elif lowername == 'nwc-eom-ccsdtq':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'energy'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'gradient'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TCE']['value'] = 'hessian'
            options ['NWCHEM']['NWCHEM_TCE_MODULE']['value'] = 'ccsdtq'

    elif lowername in dft_functionals_list(): 
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'energy'
            options ['SCF']['DFT_FUNCTIONAL']['value'] = lowername[4:].upper()
            options ['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(options) 
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'gradient'
            options ['SCF']['DFT_FUNCTIONAL']['value'] = lowername[4:].upper()
            options ['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(options)
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_DFT']['value'] = 'hessian'
            options ['SCF']['DFT_FUNCTIONAL']['value'] = lowername[4:].upper()
            options ['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(options)

    elif lowername == 'nwc-tddft':
        if dertype == 0:
            options ['NWCHEM']['NWCHEM_TASK_TDDFT']['value'] = 'energy'
        elif dertype == 1:
            options ['NWCHEM']['NWCHEM_TASK_TDDFT']['value'] = 'gradient'
        elif dertype == 2:
            options ['NWCHEM']['NWCHEM_TASK_TDDFT']['value'] = 'hessian'
 
    else:
        raise ValidationError("""Requested NWChem computational methods %s is not available.""" % (lowername))

    #Set clobbering 
    if 'NWCHEM_TASK_SCF' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_SCF']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_SCF']['superclobber'] = True

    if 'NWCHEM_TASK_DFT' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_DFT']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_DFT']['superclobber'] = True

    if 'NWCHEM_TASK_MP2' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_MP2']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_MP2']['superclobber'] = True 

    if 'NWCHEM_TASK_CCSD' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_CCSD']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_CCSD']['superclobber'] = True

    if 'NWCHEM_TASK_CCSD(T)' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_CCSD(T)']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_CCSD(T)']['superclobber'] = True 

    if 'NWCHEM_TASK_TCE' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_TCE']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_TCE']['superclobber'] = True 

    if 'NWCHEM_TASK_TDDFT' in options['NWCHEM']:
        options['NWCHEM']['NWCHEM_TASK_TDDFT']['clobber'] = True
        options['NWCHEM']['NWCHEM_TASK_TDDFT']['superclobber'] = True

    return text,options 

def harvest(p4Mol, nwout, **largs): #check orientation and scratch files
    """Parses all the pieces of output from NWChem: the stdout in
    *nwout* Scratch files are not yet considered at this moment. 

    """

    outPsivar, outMol, outGrad, version, error = harvest_output(nwout)
   # print ('outMol', outMol)
   #print ('outGrad', outGrad)
   # print ('P4Mol', p4Mol)
   #Coordinate info check

    if outMol:
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""NWChem outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
                            (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))  
            ## TEST##########
            #else: 
            #   print ( """NWChem outfile (NRE: %f) consistent with Psi4 input (NRE: %f).""" % \
                #             (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
            ################
    else:
        raise ValidationError("""No coordinate information extracted from NWChem output.""")

    return outPsivar, None, outGrad, outMol, version, error


def nwchem_list():
    """Return an array of NWChem methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-dft')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd(t)')
    val.append('nwc-eom-ccsd')
    val.append('nwc-eom-ccsdt')
    val.append('nwc-eom-ccsdtq')
    val.extend(dft_functionals_list())
    val.append('nwc-tddft')

    return val

def nwchem_gradient_list(): 
    """Return an array of NWChem methods with energies. Appended
    to procedures['gradient'].

    """
    val = []
    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-dft')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd(t)')
    val.append('nwc-eom-ccsd')
    val.append('nwc-eom-ccsdt')
    val.append('nwc-eom-ccsdtq')
    val.extend(dft_functionals_list())
    val.append('nwc-tddft')

    return val

def nwchem_hessian_list(): 
    return []


def nwchem_psivar_list():
    """Return a dict with keys of most NWChem methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.
        
    """
    VARH = {}
    VARH['nwc-scf'] = {
                            'nwc-scf' : 'HF TOTAL ENERGY'}
    VARH['nwc-hf'] = {
                             'nwc-hf' : 'HF TOTAL ENERGY'}
    VARH['nwc-mp2'] = {
                             'nwc-hf' : 'HF TOTAL ENERGY',
                            'nwc-mp2' : 'MP2 TOTAL ENERGY'}
    VARH['nwc-ccsd'] = {
                             'nwc-hf' : 'HF TOTAL ENERGY',
                            'nwc-mp2' : 'MP2 TOTAL ENERGY',
                           'nwc-ccsd' : 'CCSD TOTAL ENERGY'}
                        
    VARH['nwc-ccsdt'] = {
                              'nwc-hf' : 'HF TOTAL ENERGY',
                           'nwc-ccsdt' : 'CCSDT TOTAL ENERGY'}
    
    VARH['nwc-ccsdtq'] = {
                               'nwc-hf' : 'HF TOTAL ENERGY',
                           'nwc-ccsdtq' : 'CCSDTQ TOTAL ENERGY'}

    VARH['nwc-ccsd(t)'] = {
                             'nwc-hf' : 'HF TOTAL ENERGY',
                            'nwc-mp2' : 'MP2 TOTAL ENERGY',
                           'nwc-ccsd' : 'CCSD TOTAL ENERGY',
                        'nwc-ccsd(t)' : 'CCSD(T) TOTAL ENERGY'}

    #VARH['nwc-dft'] = {
    #                        'nwc-dfttot' : 'DFT TOTAL ENERGY'} 

    VARH['nwc-eom-ccsd'] = {
                            'nwc-hf' : 'HF TOTAL ENERGY',
                       #'nwc-dfttot' : 'DFT TOTAL ENERGY', 
                          'nwc-ccsd' : 'CCSD TOTAL ENERGY'}
    return VARH

def dft_functionals_list():
    """Return an array of available method of dft functionals through nwchem interface. 

    """
    val = []
    val.append('nwc-b3lyp')
    val.append('nwc-pbe0') 
    val.append('nwc-m05')
    val.append('nwc-m05-2x')
    val.append('nwc-ft97')
    val.append('nwc-hcth')
    val.append('nwc-hcth120')
    val.append('nwc-dldf')
    val.append('nwc-b2plyp')
    val.append('nwc-blyp')
    val.append('nwc-bp86')
    val.append('nwc-pbe')
    val.append('nwc-b97-0')
    val.append('nwc-b97-1')
    val.append('nwc-b97-2')
    val.append('nwc-hcthp14')
    val.append('nwc-hcth407p')
    val.append('nwc-b1b95')
    val.append('nwc-m06')
    val.append('nwc-m06-2x')
    val.append('nwc-m06-hf')
    val.append('nwc-m08-hx')
    val.append('nwc-m08-so')
    val.append('nwc-m11')
    val.append('nwc-m11-l')
    val.append('nwc-mpw1b95')
    val.append('nwc-mpw1k')
    val.append('nwc-mpwb1k')
    val.append('nwc-pw6b95')
    val.append('nwc-pwb6k')
    val.append('nwc-b97-gga1')
    val.append('nwc-b97')
    val.append('nwc-tpssh')
    val.append('nwc-bhandh')
    val.append('nwc-cam-b3lyp')
    val.append('nwc-bop')
#    val.append('nwc-hf')
    val.append('nwc-pw91')
    val.append('nwc-b1lyp')
    val.append('nwc-b1pw91')
    val.append('nwc-b3lyp5')
    val.append('nwc-b86bpbe')
    val.append('nwc-b97-1p')
    val.append('nwc-b97-3')
    val.append('nwc-bhandhlyp')
    val.append('nwc-bhhlyp')
    val.append('nwc-lrc-wpbe')
    val.append('nwc-lrc-wpbeh')
    val.append('nwc-mpw1pw')
    val.append('nwc-mpwlyp1m')
    val.append('nwc-mpwlyp1w')
    val.append('nwc-pbe0-13')
    val.append('nwc-pbeh')
    val.append('nwc-pbelyp1w')
    val.append('nwc-pw86pbe')
    val.append('nwc-tpsslyp1w')
    val.append('nwc-xlyp')
    val.append('nwc-wpbe')

    return val


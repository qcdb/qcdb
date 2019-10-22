import struct

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule

#   Ways of modifying a computation
#   global:     set global c-side option
#   local:      set local c-side option
#   kwarg:      set kwarg
#   i-local:    set global=local c-side option to an interface module
#   ro-def:     code uses default entirely specified by read_options
#   module-def: code uses default that is complex mixture of read_options settings
#   i-def:      interfaced code uses defaults not entirely expressed in read_options
#   driver-def: driver code sets complex defaults
#
#   Pure psi4 operation
#   kwarg ~= local > global > driver-def > module-def > ro-def
#
#   Interfaced psi4 operation
#   kwarg ~= i-local > local > global > driver-def > i-def

#   P4 infrastructure replacing interfaced infrastructure (mol, basis, mem) where unavoidable overlap in how things are specified (mult in mol{} vs keyword) is treated as a clobber & complain if conflict VS P4 infrastructure as an aliased/convenient leak into interfaced infrastructure (psi) and is strictly no clobber or complain.

# Philosophy break:
#   Specification options
#   Massaging options

#   * No program's defaults should be tampered with w/o provocation

#   want all defaults applied to all programs, so p4 scf_conv is 5 and c4 scf_conv is 5
#   want separate regimes, so conv 6 covers all the p4 parts and cfour_conv = 8 covers the c4 parts
#   want mixture, so basis gets applied to c4 but others don't
#   first case, when options specified explicitly

#   [scf][d_convergence]    [cfour][cfour_scf_conv]     what happens?
#   8 from opt()            7 by default
#   6 from set {...}        7 by default                6 (guideline that psi4 format converts when clear)
#   8 from opt()            5 from set {...}            5 (local trumps)
#   6 from set {...}        5 from set {...}            5 (local trumps)
#
#   energy(name)            [cfour][cfour_calc_level]
#   c4-scf                  SCF by default
#   c4-scf                  CCSD from set {...}


def harvest_zmat(zmat: str) -> Molecule:
    """Parses the contents of the Cfour ZMAT file into array and
    coordinate information. The coordinate info is converted into a
    rather dinky Molecule (no fragment, but does read charge, mult,
    unit). Return qcdb.Molecule. Written for findif zmat* where
    geometry always Cartesian and Bohr.

    """
    zmat = zmat.splitlines()[1:]  # skip comment line
    Nat = 0
    readCoord = True
    isBohr = ''
    charge = 0
    mult = 1
    molxyz = ''
    for line in zmat:
        if line.strip() == '':
            readCoord = False
        elif readCoord:
            molxyz += line + '\n'
            Nat += 1
        else:
            if line.find('CHARGE') > -1:
                idx = line.find('CHARGE')
                charge = line[idx + 7:]
                idxc = charge.find(',')
                if idxc > -1:
                    charge = charge[:idxc]
                charge = int(charge)
            if line.find('MULTIPLICITY') > -1:
                idx = line.find('MULTIPLICITY')
                mult = line[idx + 13:]
                idxc = mult.find(',')
                if idxc > -1:
                    mult = mult[:idxc]
                mult = int(mult)
            if line.find('UNITS=BOHR') > -1:
                isBohr = ' bohr'

    molxyz = f'{Nat}{isBohr}\n{charge} {mult}\n' + molxyz
    mol = Molecule(validate=False,
                   **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz,
                                                                       dtype='xyz+',
                                                                       fix_com=True,
                                                                       fix_orientation=True)["qm"],
                                             dtype=2))

    return mol


def format_fjobarc(energy, amap, elem, coordinates, gradient, dipole):
    """Takes the key results from a gradient computation (*energy*,
    element Z list *elem*, *coordinates*, *gradient*,
    *dipole*, and atom ordering *amap*) and writes a string *fja*
    that exactly mimics the contents of a Cfour FJOBARC file.

    """
    fja = 'TOTENERG\n'
    fja += '%15d%15d\n' % (struct.unpack("ii", struct.pack("d", energy)))
    fja += 'COORD\n'
    Nat = len(coordinates)
    flatcoord = []
    for at in range(Nat):
        for xyz in range(3):
            flatcoord.append(coordinates[amap[at]][xyz])
    for idx in range(len(flatcoord)):
        if abs(flatcoord[idx]) < 1.0E-14:  # TODO
            flatcoord[idx] = 0.0
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("d", flatcoord[idx])))
        if idx % 2 == 1:
            fja += '\n'
    if len(flatcoord) % 2 == 1:
        fja += '\n'
    fja += 'MAP2ZMAT\n'
    for idx in range(Nat):
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("l", amap[idx] + 1)))
        if idx % 2 == 1:
            fja += '\n'
    if Nat % 2 == 1:
        fja += '\n'
    fja += 'GRD FILE\n'
    fja += '%5d%20.10f\n' % (Nat, 0.0)
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], coordinates[at][0], coordinates[at][1],
                                                   coordinates[at][2])
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], gradient[at][0], gradient[at][1], gradient[at][2])
    fja += 'DIPOL FILE\n'
    fja += '%20.10f%20.10f%20.10f\n' % (dipole[0], dipole[1], dipole[2])

    return fja


def backtransform(chgeMol: 'Molecule', permMol: 'Molecule', chgeGrad=None, chgeDip=None):
    """Return `chgeMol` and `chgeGrd` to the native Cfour orientation embodied by `permMol`. Currently for vpt2.

<         p4c4 = OrientMols(p4Mol, outMol)
<         oriCoord = p4c4.transform_coordinates2(outMol)
---
>         # TODO watch out - haven't seen atom_map=False yet  # May 2019 we have now
>         rmsd, mill, amol = outMol.B787(p4Mol, atoms_map=True, mols_align=True, verbose=0)
>         oriCoord = mill.align_coordinates(outMol.geometry(np_out=True))


       rmsd, mill, amol = grdMol.B787(p4Mol, atoms_map=False, mols_align=True, verbose=0)

        oriCoord = mill.align_coordinates(grdMol.geometry(np_out=True))
        oriGrad = mill.align_gradient(np.array(grdGrad))
        if dipolDip is None:
            oriDip = None
        else:
            oriDip = mill.align_vector(np.array(dipolDip))

        if fcmHess is None:
            oriHess = None
        else:
            oriHess = mill.align_hessian(np.array(fcmHess))


    """
    # Set up array reorientation object -- opposite than usual
    amol, data = chgeMol.align(permMol, atoms_map=False, mols_align=True, verbose=0)
    mill = data['mill']

    oriCoord = mill.align_coordinates(chgeMol.geometry)
    oriElem = mill.align_atoms(np.array(chgeMol.atomic_numbers))
    oriElemMap = mill.atommap

    oriGrad = None if chgeGrad is None else mill.align_gradient(np.array(chgeGrad))
    oriDip = None if chgeDip is None else mill.align_vector(np.array(chgeDip))

    if chgeGrad is not None and chgeDip is not None:
        return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
    else:
        return oriElemMap, oriElem, oriCoord


#def backtransform_grad(p4Mol, c4Mol, p4Grd, p4Dip):
#    """Here, p4Mol and p4Grd need to be turned into the native Cfour
#    orientation embodied by c4Mol. Currently for vpt2.
#
#    """
#    # Set up array reorientation object
#    p4c4 = OrientMols(c4Mol, p4Mol)  # opposite than usual
#    oriCoord = p4c4.transform_coordinates2(p4Mol)
#    oriGrad = p4c4.transform_gradient(p4Grd)
#    p4Elem = []
#    for at in range(p4Mol.natom()):
#        p4Elem.append(p4Mol.Z(at))
#    oriElem = p4c4.transform_elementlist(p4Elem)
#    oriElemMap = p4c4.Catommap
#    oriDip = p4c4.transform_vector(p4Dip)
#
#    #print p4c4
#    #print '    <<<   Input C4 Mol   >>>'
#    #c4Mol.print_out()
#    #print '    <<<   Input P4 Mol   >>>'
#    #p4Mol.print_out()
#    #print '    <<<   Input P4 Grad   >>>'
#    #if p4Grd is not None:
#    #    for item in p4Grd:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Coord   >>>'
#    #if oriCoord is not None:
#    #    for item in oriCoord:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Elem   >>>'
#    #if oriElem is not None:
#    #    for item in oriElem :
#    #        print('       %16.8f' % (item))
#    #print '    <<<   Rotated P4 Dip  >>>'
#    #if oriDip is not None:
#    #    print('       %16.8f %16.8f %16.8f' % (oriDip[0], oriDip[1], oriDip[2]))
#    #print '    <<<   Rotated P4 Grad   >>>'
#    #if oriGrad is not None:
#    #    for item in oriGrad:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#
#    return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
#    #return oriElem, oriCoord, oriGrad, oriElemMap, oriDip


def jajo2mol(jajodic):
    """Returns a Molecule from entries in dictionary *jajodic* extracted
    from JAINDX and JOBARC.

    """
    zmap = jajodic[b'MAP2ZMAT']
    elem = jajodic[b'ATOMCHRG']
    coord = jajodic[b'COORD   ']
    Nat = len(elem)

    molxyz = '%d bohr\n\n' % (Nat)
    # TODO chgmult, though not really necessary for reorientation
    for at in range(Nat):
        posn = zmap[at] - 1
        el = 'GH' if elem[posn] == 0 else qcel.periodictable.to_E(elem[posn])
        posn *= 3
        molxyz += '%s %21.15f %21.15f %21.15f\n' % (el, coord[posn], coord[posn + 1], coord[posn + 2])
    mol = Molecule(validate=False,
                   **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz,
                                                                       dtype='xyz+',
                                                                       fix_com=True,
                                                                       fix_orientation=True)["qm"],
                                             dtype=2))

    return mol

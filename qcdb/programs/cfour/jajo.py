import struct
from typing import Dict, List


def getrec(reclabelarray: List[bytes], bjobarc: bytes, bjaindx: bytes, *, verbose: bool = False) -> Dict[bytes, bytes]:
    """Reads binary files JOBARC and JAINDX and returns contents
    of each record in *reclabelarray*.

    """
    knownlabels = {
        b"AU_LENGT": "DOUBLE",
        b"CHARGE_E": "DOUBLE",
        b"AMU     ": "DOUBLE",
        b"NUC_MAGN": "DOUBLE",
        b"MASS_ELE": "DOUBLE",
        b"MASS_PRO": "DOUBLE",
        b"HBAR    ": "DOUBLE",
        b"AU_MASSP": "DOUBLE",
        b"SP_LIGHT": "DOUBLE",
        b"AU_EV   ": "DOUBLE",
        b"AVOGADRO": "DOUBLE",
        b"AU_ENERG": "DOUBLE",
        b"AU_CM-1 ": "DOUBLE",
        b"CM-1_KCA": "DOUBLE",
        b"CM-1_KJ ": "DOUBLE",
        b"AU_DIPOL": "DOUBLE",
        b"AU_VELOC": "DOUBLE",
        b"AU_TIME ": "DOUBLE",
        b"EL_GFACT": "DOUBLE",
        b"EA_IRREP": "INTEGER",
        b"UHFRHF  ": "INTEGER",
        b"IFLAGS  ": "INTEGER",
        b"IFLAGS2 ": "INTEGER",
        b"OCCUPYA ": "INTEGER",
        b"NUMDROPA": "INTEGER",
        b"JODAFLAG": "INTEGER",
        b"TITLE   ": "CHARACTER",
        b"NCNSTRNT": "INTEGER",
        b"ICNSTRNT": "INTEGER",
        b"VCNSTRNT": "DOUBLE",
        b"NMPROTON": "INTEGER",
        b"NREALATM": "INTEGER",
        b"COORDINT": "DOUBLE",
        b"VARNAINT": "DOUBLE",
        b"COORD000": "DOUBLE",
        b"ROTCONST": "DOUBLE",
        b"ORIENT2 ": "DOUBLE",  # input orientation into interial frame
        b"LINEAR  ": "INTEGER",
        b"NATOMS  ": "INTEGER",
        b"COORD   ": "DOUBLE",
        b"ORIENTMT": "DOUBLE",  # input orientation from ZMAT (mostly useful for Cartesians) to Cfour standard orientation
        b"ATOMMASS": "DOUBLE",
        b"ORIENT3 ": "DOUBLE",
        b"FULLPTGP": "CHARACTER",
        b"FULLORDR": "INTEGER",
        b"FULLNIRR": "INTEGER",
        b"FULLNORB": "INTEGER",
        b"FULLSYOP": "DOUBLE",
        b"FULLPERM": "INTEGER",
        b"FULLMEMB": "INTEGER",
        b"FULLPOPV": "INTEGER",
        b"FULLCLSS": "INTEGER",
        b"FULLSTGP": "CHARACTER",
        b"ZMAT2MOL": "INTEGER",
        b"COMPPTGP": "CHARACTER",
        b"COMPORDR": "INTEGER",
        b"COMPNIRR": "INTEGER",
        b"COMPNORB": "INTEGER",
        b"COMPSYOP": "DOUBLE",
        b"COMPPERM": "INTEGER",
        b"COMPMEMB": "INTEGER",
        b"COMPPOPV": "INTEGER",
        b"COMPCLSS": "INTEGER",
        b"COMPSTGP": "CHARACTER",
        b"BMATRIX ": "DOUBLE",
        b"NUCREP  ": "DOUBLE",
        b"TIEDCORD": "INTEGER",
        b"MPVMZMAT": "INTEGER",
        b"ATOMCHRG": "INTEGER",
        b"NTOTSHEL": "INTEGER",
        b"NTOTPRIM": "INTEGER",
        b"BASISEXP": "DOUBLE",
        b"BASISCNT": "DOUBLE",
        b"SHELLSIZ": "INTEGER",
        b"SHELLPRM": "INTEGER",
        b"SHELLANG": "INTEGER",
        b"SHELLLOC": "INTEGER",
        b"SHOFFSET": "INTEGER",
        b"SHELLORB": "INTEGER",
        b"PROFFSET": "INTEGER",
        b"PRIMORBT": "INTEGER",
        b"FULSHLNM": "INTEGER",
        b"FULSHLTP": "INTEGER",
        b"FULSHLSZ": "INTEGER",
        b"FULSHLAT": "INTEGER",
        b"JODAOUT ": "INTEGER",
        b"NUMIIII ": "INTEGER",
        b"NUMIJIJ ": "INTEGER",
        b"NUMIIJJ ": "INTEGER",
        b"NUMIJKL ": "INTEGER",
        b"NBASTOT ": "INTEGER",
        b"NAOBASFN": "INTEGER",
        b"NUMBASIR": "INTEGER",
        b"FAOBASIR": "DOUBLE",
        b"AO2SO   ": "DOUBLE",
        b"FULLSOAO": "DOUBLE",
        b"FULLAOSO": "DOUBLE",
        b"AO2SOINV": "DOUBLE",
        b"CART3CMP": "DOUBLE",
        b"CART2CMP": "DOUBLE",
        b"CMP3CART": "DOUBLE",
        b"CMP2CART": "DOUBLE",
        b"ANGMOMBF": "INTEGER",
        b"NBASATOM": "INTEGER",
        b"NAOBFORB": "INTEGER",
        b"MAP2ZMAT": "INTEGER",
        b"CENTERBF": "INTEGER",
        b"CNTERBF0": "INTEGER",
        b"ANMOMBF0": "INTEGER",
        b"CMP2ZMAT": "DOUBLE",
        b"ZMAT2CMP": "DOUBLE",
        b"OVERLAP ": "DOUBLE",
        b"ONEHAMIL": "DOUBLE",
        b"AOOVRLAP": "DOUBLE",
        b"SHALFMAT": "DOUBLE",
        b"SCFEVCA0": "DOUBLE",
        b"RPPBMAT ": "DOUBLE",
        b"OCCUPYA0": "INTEGER",
        b"SYMPOPOA": "INTEGER",
        b"SYMPOPVA": "INTEGER",
        b"SCFEVLA0": "DOUBLE",
        b"SCFDENSA": "DOUBLE",
        b"FOCKA   ": "DOUBLE",
        b"SMHALF  ": "DOUBLE",
        b"EVECOAOA": "DOUBLE",
        b"ONEHMOA ": "DOUBLE",
        b"NOCCORB ": "INTEGER",
        b"NVRTORB ": "INTEGER",
        b"SCFENEG ": "DOUBLE",
        b"TOTENERG": "DOUBLE",
        b"IRREPALP": "INTEGER",
        b"OMEGA_A ": "DOUBLE",
        b"EVECAOXA": "DOUBLE",
        b"EVALORDR": "DOUBLE",
        b"EVECAO_A": "DOUBLE",
        b"EVCSYMAF": "CHARACTER",
        b"EVCSYMAC": "CHARACTER",
        b"TESTVECT": "DOUBLE",
        b"MODROPA ": "INTEGER",
        b"VRHARMON": "DOUBLE",
        b"NEWRECRD": "INTEGER",
        b"VRCORIOL": "DOUBLE",
        b"VRQUADRA": "DOUBLE",
        b"VRANHARM": "DOUBLE",
        b"REFINERT": "DOUBLE",
        b"DIDQ    ": "DOUBLE",
        b"REFCOORD": "DOUBLE",
        b"REFDIPOL": "DOUBLE",
        b"REFGRADI": "DOUBLE",
        b"REFDIPDR": "DOUBLE",
        b"REFNORMC": "DOUBLE",
        b"REFD2EZ ": "DOUBLE",
        b"REFFREQS": "DOUBLE",
        b"REFORIEN": "DOUBLE",
        b"NUSECORD": "INTEGER",
        b"NZMATANH": "INTEGER",
        b"ISELECTQ": "INTEGER",
        b"NEXTGEOM": "DOUBLE",
        b"NEXTGEO1": "DOUBLE",
        b"FCMDISPL": "DOUBLE",
        b"GRDDISPL": "DOUBLE",
        b"DPMDISPL": "DOUBLE",
        b"DIPDISPL": "DOUBLE",
        b"NMRDISPL": "DOUBLE",
        b"SRTDISPL": "DOUBLE",
        b"CHIDISPL": "DOUBLE",
        b"POLDISPL": "DOUBLE",
        b"EFGDISPL": "DOUBLE",
        b"THEDISPL": "DOUBLE",
        b"JFCDISPL": "DOUBLE",
        b"JSDDISPL": "DOUBLE",
        b"JSODISPL": "DOUBLE",
        b"JDSODISP": "DOUBLE",
        b"CUBCOUNT": "INTEGER",
        b"FCMMAPER": "DOUBLE",
        b"QPLSMINS": "INTEGER",
        b"CUBCOORD": "INTEGER",
        b"PASS1   ": "INTEGER",
        b"REFFORDR": "INTEGER",
        b"REFFSYOP": "DOUBLE",
        b"REFFPERM": "INTEGER",
        b"REFNUMIC": "INTEGER",
        b"REFAMAT ": "DOUBLE",
        b"REFTTEN ": "DOUBLE",
        b"REFLINER": "INTEGER",
        b"DIPOLMOM": "DOUBLE",
        b"POLARTEN": "DOUBLE",
        b"CHITENSO": "DOUBLE",
        b"EFGTENSO": "DOUBLE",
        b"IRREPPOP": "INTEGER",
        b"REORDERA": "INTEGER",
        b"IRREPBET": "INTEGER",
        b"SCFEVLB0": "DOUBLE",
        b"SCFEVCB0": "DOUBLE",
        b"IRREPCOU": "INTEGER",
        b"IDROPA  ": "INTEGER",
        b"OCCSCF  ": "INTEGER",
        b"VRTSCF  ": "INTEGER",
        b"SCFEVECA": "DOUBLE",
        b"NCOMPA  ": "INTEGER",
        b"NBASCOMP": "INTEGER",
        b"SCFEVALA": "DOUBLE",
        b"SCFEVALB": "DOUBLE",
        b"SVAVA0  ": "INTEGER",
        b"SVAVA0X ": "INTEGER",
        b"SVAVA0I ": "INTEGER",
        b"SVBVB0  ": "INTEGER",
        b"SVBVB0X ": "INTEGER",
        b"SVBVB0I ": "INTEGER",
        b"SOAOA0  ": "INTEGER",
        b"SOAOA0X ": "INTEGER",
        b"SOAOA0I ": "INTEGER",
        b"SOBOB0  ": "INTEGER",
        b"SOBOB0X ": "INTEGER",
        b"SOBOB0I ": "INTEGER",
        b"SVAVA1  ": "INTEGER",
        b"SVAVA1X ": "INTEGER",
        b"SVAVA1I ": "INTEGER",
        b"SVBVB1  ": "INTEGER",
        b"SVBVB1X ": "INTEGER",
        b"SVBVB1I ": "INTEGER",
        b"SOAOA1  ": "INTEGER",
        b"SOAOA1X ": "INTEGER",
        b"SOAOA1I ": "INTEGER",
        b"SOBOB1  ": "INTEGER",
        b"SOBOB1X ": "INTEGER",
        b"SOBOB1I ": "INTEGER",
        b"SVAOA2  ": "INTEGER",
        b"SVAOA2X ": "INTEGER",
        b"SVAOA2I ": "INTEGER",
        b"SVBOB2  ": "INTEGER",
        b"SVBOB2X ": "INTEGER",
        b"SVBOB2I ": "INTEGER",
        b"SOBVA2  ": "INTEGER",
        b"SOBVA2X ": "INTEGER",
        b"SOBVA2I ": "INTEGER",
        b"SVBOA2  ": "INTEGER",
        b"SVBOA2X ": "INTEGER",
        b"SVBOA2I ": "INTEGER",
        b"SVAVB2  ": "INTEGER",
        b"SVAVB2X ": "INTEGER",
        b"SVAVB2I ": "INTEGER",
        b"SOAOB2  ": "INTEGER",
        b"SOAOB2X ": "INTEGER",
        b"SOAOB2I ": "INTEGER",
        b"SOAVA2  ": "INTEGER",
        b"SOAVA2X ": "INTEGER",
        b"SOAVA2I ": "INTEGER",
        b"SOBVB2  ": "INTEGER",
        b"SOBVB2X ": "INTEGER",
        b"SOBVB2I ": "INTEGER",
        b"SOAVB2  ": "INTEGER",
        b"SOAVB2X ": "INTEGER",
        b"SOAVB2I ": "INTEGER",
        b"SVAVA2  ": "INTEGER",
        b"SVAVA2X ": "INTEGER",
        b"SVAVA2I ": "INTEGER",
        b"SVBVB2  ": "INTEGER",
        b"SVBVB2X ": "INTEGER",
        b"SVBVB2I ": "INTEGER",
        b"SOAOA2  ": "INTEGER",
        b"SOAOA2X ": "INTEGER",
        b"SOAOA2I ": "INTEGER",
        b"SOBOB2  ": "INTEGER",
        b"SOBOB2X ": "INTEGER",
        b"SOBOB2I ": "INTEGER",
        b"SYMPOPOB": "INTEGER",
        b"SYMPOPVB": "INTEGER",
        b"T2NORM  ": "DOUBLE",
        b"MOIOVEC ": "INTEGER",
        b"MOIOWRD ": "INTEGER",
        b"MOIOSIZ ": "INTEGER",
        b"MOIODIS ": "INTEGER",
        b"MOIOFIL ": "INTEGER",
        b"ISYMTYP ": "INTEGER",
        b"TOTRECMO": "INTEGER",
        b"TOTWRDMO": "INTEGER",
        b"RELDENSA": "DOUBLE",
        b"IINTERMA": "DOUBLE",
        b"OCCNUM_A": "DOUBLE",
        b"SCRATCH ": "DOUBLE",
        b"SETUP2  ": "INTEGER",
        b"MOLHES2 ": "INTEGER",
        b"GRAD2   ": "INTEGER",
        b"COORDMAS": "INTEGER",
        b"NUCMULT ": "INTEGER",
        b"SYMCOORD": "DOUBLE",
        b"SYMCOOR2": "DOUBLE",
        b"SYMCOOR3": "DOUBLE",
        b"SYMMLENG": "INTEGER",
        b"SKIP    ": "INTEGER",
        b"NSYMPERT": "INTEGER",
        b"NPERTB  ": "INTEGER",
        b"TRANSINV": "INTEGER",
        b"IBADNUMB": "INTEGER",
        b"IBADINDX": "INTEGER",
        b"IBADIRRP": "INTEGER",
        b"IBADPERT": "INTEGER",
        b"IBADSPIN": "INTEGER",
        b"TREATPER": "INTEGER",
        b"MAXAODSZ": "INTEGER",
        b"PERTINFO": "INTEGER",
        b"GRADIENT": "DOUBLE",
        b"HESSIANM": "DOUBLE",
        b"GRDZORDR": "DOUBLE",
        b"D2EZORDR": "DOUBLE",
        b"REALCORD": "DOUBLE",
        b"DUMSTRIP": "INTEGER",
        b"BMATRIXC": "DOUBLE",
        b"REALATOM": "INTEGER",
        b"NORMCORD": "DOUBLE",
        b"DIPDERIV": "DOUBLE",
        b"I4CDCALC": "DOUBLE",
        b"FREQUENC": "DOUBLE",
        b"RATMMASS": "DOUBLE",
        b"RATMPOSN": "INTEGER",
        b"DEGENERT": "INTEGER",
        b"REFSHILD": "DOUBLE",
        b"CORIZETA": "DOUBLE",
        b"NMPOINTX": "INTEGER",
        b"REFD3EDX": "DOUBLE",
        b"BPPTOB  ": "DOUBLE",
        b"BPTOB   ": "DOUBLE",
        b"BSRTOB  ": "DOUBLE",
        b"BARTOB  ": "DOUBLE",
        b"VRTOTAL ": "DOUBLE",
        b"D2DIPOLE": "DOUBLE",
        b"D3DIPOLE": "DOUBLE",
        b"D1DIPOLE": "DOUBLE",
        b"REFNORM2": "DOUBLE",
        b"NUSECOR2": "INTEGER",
        b"FCMDISP2": "DOUBLE",
        b"RGTDISPL": "DOUBLE",
        b"CUBCOOR1": "INTEGER",
        b"CUBCOOR2": "INTEGER",
        b"REFFPEM2": "INTEGER",
        b"RGTTENSO": "DOUBLE",
        b"REFFPER2": "INTEGER",
        b"REFD4EDX": "DOUBLE",
        b"ZPE_ANHA": "DOUBLE",
        b"OPENSLOT": "INTEGER",
        b"BOLTZMAN": "DOUBLE",
        b"MRCCOCC ": "INTEGER",
        b"ABELPTGP": "CHARACTER",
        b"ABELORDR": "INTEGER",
        b"ABELNIRR": "INTEGER",
        b"ABELNORB": "INTEGER",
        b"ABELSYOP": "DOUBLE",
        b"ABELPERM": "INTEGER",
        b"ABELMEMB": "INTEGER",
        b"ABELPOPV": "INTEGER",
        b"ABELCLSS": "INTEGER",
        b"ABELSTGP": "CHARACTER",
        b"REALCHRG": "INTEGER",  # atom/mol? charge taking into acct edp
        b"NSOSCF  ": "INTEGER",  # whether is spin orbital calc?
        b"SCFVCFLA": "DOUBLE",  # scf vector expanded from sph to cart basis for symm anal - determin orb sym
        b"EFG_SYM1": "INTEGER",  # symmetry property of components of electric field gradient  integrals
        b"EFG_SYM2": "INTEGER",  # symm prop of comp of EFG
        b"DCTDISPL": "DOUBLE",
        b"DANGERUS": "INTEGER",  # ?
        b"FULLCHAR": "CHARACTER",  # ?
        b"FULLDEGN": "CHARACTER",  # ?
        b"FULLLABL": "CHARACTER",  # ?
        b"FULLNIRX": "CHARACTER",  # ?
        b"COMPCHAR": "CHARACTER",  # ?
        b"COMPDEGN": "CHARACTER",  # ?
        b"COMPLABL": "CHARACTER",  # ?
        b"COMPNIRX": "CHARACTER",  # ?
        b"ROTVECX ": "CHARACTER",  # ?
        b"ROTVECY ": "CHARACTER",  # ?
        b"ROTVECZ ": "CHARACTER",  # ?
        b"COMPNSYQ": "CHARACTER",  # ?
        b"COMPSYQT": "CHARACTER",  # ?
        b"COMPSYMQ": "CHARACTER",  # ?
        b"TRAVECX ": "CHARACTER",  # ?
        b"TRAVECY ": "CHARACTER",  # ?
        b"TRAVECZ ": "CHARACTER",  # ?
        b"NVIBSYM ": "CHARACTER",  # ?
        b"NUMVIBRT": "CHARACTER",  # ?
        b"SBGRPSYM": "CHARACTER",  # ?
        b"ORDERREF": "CHARACTER",  # ?
        b"OPERSREF": "CHARACTER",  # ?
        b"NVIBSYMF": "CHARACTER",  # ?
        b"FULLNSYQ": "CHARACTER",  # ?
        b"FULLSYQT": "CHARACTER",  # ?
        b"FULLSYMQ": "CHARACTER",  # ?
        b"INVPSMAT": "CHARACTER",  # ?
        b"FDCOORDS": "CHARACTER",  # ?
        b"FDCALCTP": "CHARACTER",  # ?
        b"NUMPOINT": "CHARACTER",  # ?
        b"NPTIRREP": "CHARACTER",  # ?
        b"GRDPOINT": "CHARACTER",  # ?
        b"DIPPOINT": "CHARACTER",  # ?
        b"ENGPOINT": "CHARACTER",  # ?
        b"PASS1FIN": "CHARACTER",  # ?
        b"REFENERG": "CHARACTER",  # ?
        b"NEXTCALC": "CHARACTER",  # ?
        b"PRINSPIN": "CHARACTER",  # ?
        b"PRINFROM": "CHARACTER",  # ?
        b"PRININTO": "CHARACTER",  # ?
        b"NEXTGEOF": "CHARACTER",  # ?
        b"ZPE_HARM": "DOUBLE",  # ?
        b"NDROPPED": "INTEGER",
        b"REFCPTGP": "INTEGER",  # ?
        b"REFFPTGP": "INTEGER",  # ?
    }

    # with open('JAINDX', mode='rb') as file:  # b is important -> binary
    fileContent = bjaindx  # file.read()
    fileLength = len(fileContent)

    if fileLength == 16012:
        srcints = 4
        srcrecs = 4
    elif fileLength == 16020:
        srcints = 4
        srcrecs = 8
    elif fileLength == 24016:
        srcints = 8
        srcrecs = 4
    elif fileLength == 24024:
        srcints = 8
        srcrecs = 8

    # fixed number of slots for options
    nopt = 1000

    type2len = {
        "DOUBLE": 8,
        "INTEGER": srcints,
        "CHARACTER": 1,
    }

    intlen2format = {
        4: "i",
        8: "l",
    }

    type2format = {
        "DOUBLE": "d",
        "INTEGER": intlen2format[type2len["INTEGER"]],
        "CHARACTER": "c",
    }

    if verbose:
        print("\n<<<  JAINDX  >>>\n")

    posf = srcrecs
    istr = intlen2format[srcrecs]
    jastart = struct.unpack(istr, fileContent[:posf])
    if verbose:
        print("%10s%10d%10d" % ("start", 0, posf))

    poss = posf
    posf = poss + 8 * nopt
    istr = "8s" * nopt
    jaindx = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print("%10s%10d%10d" % ("jaindx", poss, posf))

    poss = posf
    posf = poss + srcints * nopt
    istr = intlen2format[srcints] * nopt
    jaindx2 = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print("%10s%10d%10d" % ("jaindx2", poss, posf))

    poss = posf
    posf = poss + srcints * nopt
    istr = intlen2format[srcints] * nopt
    jaindx3 = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print("%10s%10d%10d" % ("jaindx3", poss, posf))

    poss = posf
    posf = poss + srcints
    istr = intlen2format[srcints]
    jamid = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print("%10s%10d%10d" % ("mid", poss, posf))

    poss = posf
    posf = poss + srcrecs
    istr = intlen2format[srcrecs]
    jaend = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print("%10s%10d%10d" % ("end", poss, posf))

    nrecs = jaindx.index(b"OPENSLOT")  # number of active records

    if verbose:
        print("\n")
        print("%20s%10d" % ("File Length:", fileLength))
        print("%20s%10d" % ("srcints Int Length:", srcints))
        print("%20s%10d" % ("srcrecs Int Length:", srcrecs))
        print("%20s%10d" % ("First Rec:", jastart[0]))
        print("%20s%10d" % ("Second Rec:", jamid[0]))
        print("%20s%10d" % ("Last Rec:", jaend[0]))
        print("%20s%10d" % ("Full Records:", nrecs))
        print("\n")

        print("\n<<<  JOBARC  >>>\n")

    # with open('JOBARC', mode='rb') as file:  # b is important -> binary
    fileContent = bjobarc  # file.read()

    returnRecords = {}
    poss = 0
    for item in range(nrecs):
        posf = poss + type2len[knownlabels[jaindx[item]]] * jaindx3[item]
        istr = type2format[knownlabels[jaindx[item]]] * jaindx3[item]
        if knownlabels[jaindx[item]] == "CHARACTER":
            bound = type2len[knownlabels[jaindx[item]]] * jaindx3[item] * 8
            posf = poss + bound
            istr = str(bound) + "s"
        jobarc = struct.unpack(istr, fileContent[poss:posf])

        if verbose:
            # print item, istr, poss, posf, '\t', jaindx[item], jaindx2[item], jaindx3[item], jobarc
            if jaindx3[item] < 120:
                print(jaindx[item], jaindx2[item], jaindx3[item], jobarc)

        poss = posf
        if jaindx[item] in reclabelarray:
            returnRecords[jaindx[item]] = jobarc

    return returnRecords


# if __name__ == "__main__":
#    want = ['NATOMS  ', 'AU_LENGT', 'COORD   ', 'HBAR    ', 'ATOMCHRG']
##    got = get_jajo_record(want)
#    got = getrec(want)
#    for item in got.keys():
#        print item, got[item]

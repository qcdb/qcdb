import sys
import uuid


# @moptions.register_opts(pe.nu_options)
def muster_inherited_keywords(ropts: "Keywords", mcfg: "ModeConfig", verbose: int = 1) -> None:
    accession = sys._getframe().f_code.co_name + "_" + str(uuid.uuid4())
    kwgs = {"accession": accession, "verbose": verbose}
    ropts.scroll["QCDB"]["TRANSLATE_QCDB"].value

    #    # qcdb/memory [B] --> cfour/memory_size [MB]
    #    qopt = ropts.scroll['QCDB']['MEMORY']
    #    if do_translate or qopt.is_required():
    #        mem = int(0.000001 * qopt.value)
    #        print('\n\nMEMORY', mem, '\n\n')
    #        ropts.suggest('CFOUR', 'MEMORY_SIZE', mem, **kwgs)
    #        ropts.suggest('CFOUR', 'MEM_UNIT', 'MB', **kwgs)
    #
    #    # qcdb/puream --> cfour/spherical
    #    ropts.suggest('CFOUR', 'SPHERICAL', ropts.scroll['QCDB']['PUREAM'].value, **kwgs)
    #
    #    # qcdb/reference --> cfour/reference
    #    # TODO ref or scf__ref?
    #    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    #    if qref in ['RHF', 'UHF', 'ROHF']:
    #    #ref = {'RHF': 'RHF',
    #    #       'UHF': 'UHF',
    #    #       'ROHF': 'ROHF'}[ropts.scroll['QCDB']['REFERENCE'].value]
    #        ropts.suggest('CFOUR', 'REFERENCE', qref, **kwgs)

    # qcdb/mp2__mp2_type --> psi4/mp2_type
    qopt = ropts.scroll["QCDB"]["MP2__MP2_TYPE"]
    if qopt.disputed():
        ropts.suggest("PSI4", "MP2_TYPE", qopt.value, **kwgs)

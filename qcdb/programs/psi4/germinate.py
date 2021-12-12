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

    qopt = ropts.scroll["QCDB"]["BASIS"]
    ropts.require("PSI4", "BASIS", qopt.value, **kwgs)  # require?

    # qcdb/reference --> psi4/scf__reference
    qopt = ropts.scroll["QCDB"]["REFERENCE"]
    if qopt.disputed():
        ropts.suggest("PSI4", "REFERENCE", qopt.value, **kwgs)
        ropts.suggest("PSI4", "SCF__REFERENCE", qopt.value, **kwgs)

    # qcdb/scf_type --> psi4/scf_type
    qopt = ropts.scroll["QCDB"]["SCF_TYPE"]
    if mcfg.translate_method_algorithm is True:
        val = qopt.value
    elif qopt.disputed():
        val = qopt.value2
    else:
        val = "skip"
    if val != "skip":
        if val == "CONV":
            val = "PK"

        ropts.suggest("PSI4", "SCF_TYPE", val, **kwgs)
        ropts.suggest("PSI4", "SCF__SCF_TYPE", val, **kwgs)

    # qcdb/mp2_type --> psi4/mp2_type
    qopt = ropts.scroll["QCDB"]["MP2_TYPE"]
    if mcfg.translate_method_algorithm is True:
        val = qopt.value
    elif qopt.disputed():
        val = qopt.value2
    else:
        val = "skip"
    if val != "skip":
        ropts.suggest("PSI4", "MP2_TYPE", val, **kwgs)

    # qcdb/freeze_core --> psi4/freeze_core
    qopt = ropts.scroll["QCDB"]["FREEZE_CORE"]
    val = "skip"
    if mcfg.translate_orbital_space is True:
        val = qopt.value
    elif qopt.disputed():
        val = qopt.value2

    if val != "skip":
        ropts.suggest("PSI4", "FREEZE_CORE", val, **kwgs)
        ropts.suggest("PSI4", "SCF__FREEZE_CORE", val, **kwgs)

    # qcdb/e_convergence --> psi4/e_convergence
    qopt = ropts.scroll["QCDB"]["E_CONVERGENCE"]
    if qopt.disputed():
        val = qopt.value
    else:
        val = "skip"

    if val != "skip":
        ropts.suggest("PSI4", "E_CONVERGENCE", val, **kwgs)

    # qcdb/qc_module --> psi4/qc_module
    qopt = ropts.scroll["QCDB"]["QC_MODULE"]
    if qopt.disputed():
        val = qopt.value
    else:
        val = "skip"

    if val != "skip":
        ropts.suggest("PSI4", "QC_MODULE", val, **kwgs)

"""
from https://github.com/psi4/psi4/blob/master/tests/tu1-h2o-energy/input.dat
Sample HF/cc-pVDZ H2O computation

"""
import pprint

import pytest

import qcdb

from .utils import *

tu1_scf_ene = -76.02665366
tu1_scf_ene_df = -76.0266327341
tu1_mp2_ae_ene = -76.23080846
tu1_mp2_ae_ene_df = -76.2307574465
tu1_mp2_fc_ene = -76.2284754833
tu1_mp2_fc_ene_df = -76.2284245503
conv = tu1_scf_ene
df = tu1_scf_ene_df

hf = tu1_scf_ene
dfhf = tu1_scf_ene_df
hfqcdbdef = hf
mp2 = tu1_mp2_ae_ene
dfmp2 = tu1_mp2_ae_ene_df
mp2fc = tu1_mp2_fc_ene
dfmp2fc = tu1_mp2_fc_ene_df
mp2qcdbdef = mp2  # usual
#mp2qcdbdef = mp2fc  # if change qcdb/keywords/read_options.py freeze_core default=True

_sh2o = """
  O
  H 1 0.96
  H 1 0.96 2 104.5
"""

@using("psi4")
@pytest.mark.parametrize("mode_options,keywords,xptd", [
    pytest.param({"translate_method_algorithm": True}, {}, conv),  # 0
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "pk"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "df"}, df),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "conv"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "pk"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "df"}, df),  # 5
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "pk", "scf_type": "conv"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "pk", "scf_type": "pk"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "pk", "scf_type": "df"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "df", "scf_type": "conv"}, df),
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "df", "scf_type": "pk"}, df),  # 10
    pytest.param({"translate_method_algorithm": True}, {"psi4_scf_type": "df", "scf_type": "df"}, df),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "conv", "psi4_scf_type": "pk"}, conv),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "pk", "psi4_scf_type": "pk"},  conv),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "df", "psi4_scf_type": "pk"},  conv),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "conv", "psi4_scf_type": "df"}, df),  # 15
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "pk", "psi4_scf_type": "df"},  df),
    pytest.param({"translate_method_algorithm": True}, {"scf_type": "df", "psi4_scf_type": "df"},  df),

    pytest.param({"translate_method_algorithm": False}, {}, df),                                # !=!
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "pk"}, conv),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "df"}, df),  # 20
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "conv"}, conv),            # !=!  # orig df
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "pk"}, conv),              # !=!  # orig df
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "df"}, df),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "pk", "scf_type": "conv"}, conv),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "pk", "scf_type": "pk"}, conv),  # 25
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "pk", "scf_type": "df"}, conv),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "df", "scf_type": "conv"}, df),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "df", "scf_type": "pk"}, df),
    pytest.param({"translate_method_algorithm": False}, {"psi4_scf_type": "df", "scf_type": "df"}, df),
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "conv", "psi4_scf_type": "pk"}, conv),  # 30
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "pk", "psi4_scf_type": "pk"},  conv),
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "df", "psi4_scf_type": "pk"},  conv),
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "conv", "psi4_scf_type": "df"}, df),
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "pk", "psi4_scf_type": "df"},  df),
    pytest.param({"translate_method_algorithm": False}, {"scf_type": "df", "psi4_scf_type": "df"},  df),  # 35
])
def test_mode_psi_hf_details(mode_options, keywords, xptd):
    h2o = qcdb.set_molecule(_sh2o)
    #qcdb.set_keywords({"scf_type": "pk"})
    qcdb.set_keywords(keywords)

    ene, wfn = qcdb.energy("p4-hf/cc-pVDZ", mode_options=mode_options, return_wfn=True)
    pprint.pprint(wfn, width=200)

    assert compare_values(xptd, ene, 6, "energy")


@pytest.mark.parametrize("program,mode_options,keywords,xptd", [
    pytest.param("cfour", {"translate_method_algorithm": True}, {}, hfqcdbdef, marks=using("cfour")),  # 0
    #
    #
    pytest.param("cfour", {"translate_method_algorithm": True}, {"scf_type": "conv"}, hf, marks=using("cfour")),
    pytest.param("cfour", {"translate_method_algorithm": True}, {"scf_type": "df"}, "raise", marks=using("cfour")),

    pytest.param("cfour", {"translate_method_algorithm": False}, {}, hf, marks=using("cfour")),  # 3
    #
    #
    pytest.param("cfour", {"translate_method_algorithm": False}, {"scf_type": "conv"}, hf, marks=using("cfour")),
    pytest.param("cfour", {"translate_method_algorithm": False}, {"scf_type": "df"}, "raise", marks=using("cfour")),

    pytest.param("gamess", {"translate_method_algorithm": True}, {}, hfqcdbdef, marks=using("gamess")),  # 6
    #
    #
    pytest.param("gamess", {"translate_method_algorithm": True}, {"scf_type": "conv"}, hf, marks=using("gamess")),
    pytest.param("gamess", {"translate_method_algorithm": True}, {"scf_type": "df"}, "raise", marks=using("gamess")),

    pytest.param("gamess", {"translate_method_algorithm": False}, {}, hf, marks=using("gamess")),  # 9
    #
    #
    pytest.param("gamess", {"translate_method_algorithm": False}, {"scf_type": "conv"}, hf, marks=using("gamess")),
    pytest.param("gamess", {"translate_method_algorithm": False}, {"scf_type": "df"}, "raise", marks=using("gamess")),

    pytest.param("nwchem", {"translate_method_algorithm": True}, {}, hfqcdbdef, marks=using("nwchem")),  # 12
    #
    #
    pytest.param("nwchem", {"translate_method_algorithm": True}, {"scf_type": "conv"}, hf, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_method_algorithm": True}, {"scf_type": "df"}, "raise", marks=using("nwchem")),

    pytest.param("nwchem", {"translate_method_algorithm": False}, {}, hf, marks=using("nwchem")),  # 15
    #
    #
    pytest.param("nwchem", {"translate_method_algorithm": False}, {"scf_type": "conv"}, hf, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_method_algorithm": False}, {"scf_type": "df"}, "raise", marks=using("nwchem")),

    pytest.param("psi4", {"translate_method_algorithm": True}, {}, hfqcdbdef, marks=using("psi4")),  # 18
    pytest.param("psi4", {"translate_method_algorithm": True}, {"psi4_scf_type": "pk"}, hf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": True}, {"psi4_scf_type": "df"}, dfhf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": True}, {"scf_type": "conv"}, hf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": True}, {"scf_type": "df"}, dfhf, marks=using("psi4")),

    pytest.param("psi4", {"translate_method_algorithm": False}, {}, dfhf, marks=using("psi4")),  # 23
    pytest.param("psi4", {"translate_method_algorithm": False}, {"psi4_scf_type": "pk"}, hf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": False}, {"psi4_scf_type": "df"}, dfhf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": False}, {"scf_type": "conv"}, hf, marks=using("psi4")),
    pytest.param("psi4", {"translate_method_algorithm": False}, {"scf_type": "df"}, dfhf, marks=using("psi4")),
])
def test_mode_dfconv_hf(program, mode_options, keywords, xptd):

    mode_options = {"translate_method_algorithm": True, **mode_options}
    prefix = qcdb.util.program_prefix(program)

    h2o = qcdb.set_molecule(_sh2o)
    qcdb.set_keywords(keywords)

    if xptd == "raise":
        with pytest.raises(qcdb.ValidationError) as e:
            qcdb.energy(prefix + "hf/cc-pVDZ", mode_options=mode_options, return_wfn=True)

        assert "SCF_TYPE 'DF' is not available" in str(e)

    else:
        ene, wfn = qcdb.energy(prefix + "hf/cc-pVDZ", mode_options=mode_options, return_wfn=True)
        pprint.pprint(wfn, width=200)

        assert compare_values(xptd, ene, atol=1.e-6, label="energy")


@pytest.mark.parametrize("program,mode_options,keywords,xptd", [
    pytest.param("cfour", {"translate_orbital_space": True}, {}, mp2qcdbdef, marks=using("cfour")),  # 0
    pytest.param("cfour", {"translate_orbital_space": True}, {"cfour_frozen_core": False}, mp2, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": True}, {"cfour_frozen_core": True}, mp2fc, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": True}, {"freeze_core": False}, mp2, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": True}, {"freeze_core": True}, mp2fc, marks=using("cfour")),

    pytest.param("cfour", {"translate_orbital_space": False}, {}, mp2, marks=using("cfour")),  # 5
    pytest.param("cfour", {"translate_orbital_space": False}, {"cfour_frozen_core": False}, mp2, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": False}, {"cfour_frozen_core": True}, mp2fc, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": False}, {"freeze_core": False}, mp2, marks=using("cfour")),
    pytest.param("cfour", {"translate_orbital_space": False}, {"freeze_core": True}, mp2fc, marks=using("cfour")),

    pytest.param("gamess", {"translate_orbital_space": True}, {}, mp2qcdbdef, marks=using("gamess")),  # 10
    pytest.param("gamess", {"translate_orbital_space": True}, {"gamess_mp2__nacore": 0}, mp2, marks=using("gamess")),
    pytest.param("gamess", {"translate_orbital_space": True}, {"gamess_mp2__nacore": 1}, mp2fc, marks=using("gamess")),
    pytest.param("gamess", {"translate_orbital_space": True}, {"freeze_core": False}, mp2, marks=using("gamess")),
    pytest.param("gamess", {"translate_orbital_space": True}, {"freeze_core": True}, mp2fc, marks=using("gamess")),

    pytest.param("gamess", {"translate_orbital_space": False}, {}, mp2fc, marks=using("gamess")),  # 15
    pytest.param("gamess", {"translate_orbital_space": False}, {"gamess_mp2__nacore": 0}, mp2, marks=using("gamess")),
    pytest.param("gamess", {"translate_orbital_space": False}, {"gamess_mp2__nacore": 1}, mp2fc, marks=using("gamess")),
    pytest.param("gamess", {"translate_orbital_space": False}, {"freeze_core": False}, mp2, marks=using("gamess")),  # ?
    pytest.param("gamess", {"translate_orbital_space": False}, {"freeze_core": True}, mp2fc, marks=using("gamess")),  # ?

    pytest.param("nwchem", {"translate_orbital_space": True}, {}, mp2qcdbdef, marks=using("nwchem")),  # 20
    pytest.param("nwchem", {"translate_orbital_space": True}, {"nwchem_mp2__freeze__core__atomic": False,}, mp2, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": True}, {"nwchem_mp2__freeze__core__atomic": True,}, mp2fc, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": True}, {"freeze_core": False}, mp2, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": True}, {"freeze_core": True}, mp2fc, marks=using("nwchem")),

    pytest.param("nwchem", {"translate_orbital_space": False}, {}, mp2, marks=using("nwchem")),  # 25
    pytest.param("nwchem", {"translate_orbital_space": False}, {"nwchem_mp2__freeze__core__atomic": False}, mp2, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": False}, {"nwchem_mp2__freeze__core__atomic": True}, mp2fc, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": False}, {"freeze_core": False}, mp2, marks=using("nwchem")),
    pytest.param("nwchem", {"translate_orbital_space": False}, {"freeze_core": True}, mp2fc, marks=using("nwchem")),

    pytest.param("psi4", {"translate_orbital_space": True}, {}, mp2qcdbdef, marks=using("psi4")),  # 30
    pytest.param("psi4", {"translate_orbital_space": True}, {"psi4_freeze_core": False}, mp2, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": True}, {"psi4_freeze_core": True}, mp2fc, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": True}, {"freeze_core": False}, mp2, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": True}, {"freeze_core": True}, mp2fc, marks=using("psi4")),

    pytest.param("psi4", {"translate_orbital_space": False}, {}, mp2, marks=using("psi4")),  # 35
    pytest.param("psi4", {"translate_orbital_space": False}, {"psi4_freeze_core": False}, mp2, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": False}, {"psi4_freeze_core": True}, mp2fc, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": False}, {"freeze_core": False}, mp2, marks=using("psi4")),
    pytest.param("psi4", {"translate_orbital_space": False}, {"freeze_core": True}, mp2fc, marks=using("psi4")),
])
def test_mode_fcae_mp2(program, mode_options, keywords, xptd):

    mode_options = {"translate_method_algorithm": True, **mode_options}
    prefix = qcdb.util.program_prefix(program)

    h2o = qcdb.set_molecule(_sh2o)
    qcdb.set_keywords(keywords)

    ene, wfn = qcdb.energy(prefix + "mp2/cc-pVDZ", mode_options=mode_options, return_wfn=True)
    pprint.pprint(wfn, width=200)

    assert compare_values(xptd, ene, atol=1.e-6, label="energy")


@pytest.mark.parametrize("mo,xptd", [
    pytest.param({}, {"mode": None, "implicit_program": "qcdb", "module_fallback": None, "translate_method_algorithm": None, "translate_orbital_space": None}),
    pytest.param({"translate_method_algorithm": True}, {"mode": None, "implicit_program": "qcdb", "module_fallback": None, "translate_method_algorithm": True, "translate_orbital_space": None}),
    pytest.param({"translate_method_algorithm": True, "translate_orbital_space": False}, {"mode": None, "implicit_program": "qcdb", "module_fallback": None, "translate_method_algorithm": True, "translate_orbital_space": False}),
    pytest.param({"mode": "unified"}, {"mode": "unified", "implicit_program": "qcdb", "module_fallback": True, "translate_method_algorithm": True, "translate_orbital_space": True}),
    pytest.param({"mode": "sandwich"}, {"mode": "sandwich", "implicit_program": "qcdb", "module_fallback": False, "translate_method_algorithm": False, "translate_orbital_space": False}),
    pytest.param({"mode": "unified", "translate_method_algorithm": False}, {"mode": "unified", "implicit_program": "qcdb", "module_fallback": True, "translate_method_algorithm": False, "translate_orbital_space": True}),
    pytest.param({"mode": "sandwich", "translate_orbital_space": True}, {"mode": "sandwich", "implicit_program": "qcdb", "module_fallback": False, "translate_method_algorithm": False, "translate_orbital_space": True}),
])
def test_mode_config(mo, xptd):
    mode_options = qcdb.driver.config.get_mode_config(mode_options=mo)

    print(mode_options)
    assert mode_options.dict() == xptd

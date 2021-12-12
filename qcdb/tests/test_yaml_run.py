import qcdb

from .utils import *

# perhaps from pywrap_cbs1 psi4 test case

yamlin_hf3 = """
molecule: |
  Li

driver: !!python/name:qcdb.energy
method: !!python/name:qcdb.cbs
kwargs:
  scf_wfn: pkg-hf
  scf_basis: aug-cc-pV[DTQ]Z
  scf_scheme: !!python/name:qcdb.scf_xtpl_helgaker_3
options:
  reference: rohf
  scf_type: conv  # longstanding pk
  psi4_df_scf_guess: false
"""
ans_hf3 = -7.4326961561955551


@using("psi4")
def test_hf3_a():
    subject = yamlin_hf3.replace("pkg-", "")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_hf3, qcdb.variable("current energy"), 7, "")


@using("psi4")
def test_hf3_b():
    subject = yamlin_hf3.replace("pkg", "p4")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_hf3, qcdb.variable("current energy"), 7, "")


@using("cfour")
def test_hf3_c():
    subject = yamlin_hf3.replace("pkg", "c4")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_hf3, qcdb.variable("current energy"), 7, "")


yamlin_c2d2 = """
molecule: |
  H
  H 1 R

  R = 1

driver: !!python/name:qcdb.energy
method: !!python/name:qcdb.cbs
kwargs:
  corl_wfn: pkg-mp2
  corl_basis: cc-pV[TQ]Z
  corl_scheme: !!python/name:qcdb.corl_xtpl_helgaker_2
  delta_wfn: pkg-ccsd
  delta_basis: cc-pV[DT]Z
  delta_scheme: !!python/name:qcdb.corl_xtpl_helgaker_2
options:
  reference: rhf
  psi4_guess: core
  mp2_type: conv
#  scf_type: pk
  psi4_df_scf_guess: true
"""
ans_c2d2 = -1.148287763304


@using("cfour")
def test_c2d2_a():
    subject = yamlin_c2d2.replace("pkg-", "")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_c2d2, qcdb.variable("current energy"), 7, "")


@using("cfour")
def test_c2d2_b():
    subject = yamlin_c2d2.replace("pkg", "p4")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_c2d2, qcdb.variable("current energy"), 7, "")


@using("cfour")
def test_c2d2_c():
    subject = yamlin_c2d2.replace("pkg", "c4")

    jrec = qcdb.yaml_run(subject)
    assert compare_values(ans_c2d2, qcdb.variable("current energy"), 7, "")


# TODO do mixed pkgs and assert provenance pattern

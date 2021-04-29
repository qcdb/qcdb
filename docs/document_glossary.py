import sys
from pathlib import Path

def document_glossary(outfile: str) -> None:

    path_to_qcdb = Path("../qcdb").resolve().parent
    sys.path.append(str(path_to_qcdb))
    import qcdb
    from qcdb.qcvars.glossary import qcvardefs

    rst = []
    rst.append(".. _`apdx:qcvariables_alpha`:")
    rst.append("")
    rst.append("QCVariables by Alpha")
    rst.append("====================")
    rst.append("")

    for qcvar, info in sorted(qcvardefs.items()):
        rst.append(f".. qcvar:: {qcvar}\n")

        for line in info["glossary"].split("\n"):
            if line.strip():
                rst.append(f"   {line.strip().replace('???', '   ')}")
        rst.append(f"   units: [{info['units']}]")
        if "dimension" in info:
            rst.append(f"   dimension: [{info['dimension']}]")
        rst.append("")

    with open(outfile, "w") as fp:
        fp.write("\n".join(rst))

#    print("\n".join(rst))


if __name__ == "__main__":
    document_glossary("source/autodoc_glossary_qcvars.rst")

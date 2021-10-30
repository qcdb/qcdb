import ast


def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


def fcell(winnowed):
    global notes

    if winnowed:
        if len(winnowed) > 1 and len(set([entry["status"] for entry in winnowed])) > 1:
            print("trouble:")
            for entry in winnowed:
                print("\t", entry)

        status = winnowed[0]["status"]
        note = winnowed[0]["note"]
        if status == "pass":
            cell = "Y"
            cell = u'\u2713'
#            cell = '✅'
            cell = u'\u2705'
            cell = '✅'
            cell = u'\u2713'
        elif status == "fd":
            cell = "fd"
        elif status == "error":
            cell = "---"
            cell = u'\u2717'
        elif status == "wrong":
            cell = "bad"
            cell = '💥 '
            cell = u'\u25a0'

        if note:
            notes.append(note)
            notes = unique(notes)
            idx = notes.index(note) + 10
            cell += f" [#f{idx}]_"
            
    else:
        cell = "?"
        cell = ""

    return cell


def fheader(head, width):
    trans = {
        # eghs
        "energy": [":py:func:`~qcdb.energy()`", "Energy", "E"],
        "gradient": [":py:func:`~qcdb.gradient()`", "Gradient", "G"],
        "hessian": [":py:func:`~qcdb.hessian()`", "Hessian", "H"],
        # fcae
        "ae": ["All-Electron", "AE"],
        "fc": ["Frozen-Core", "FC"],
        # refs
        "rhf": ["RHF"],
        "uhf": ["UHF"],
        "rohf": ["ROHF"],
    }
    if head in trans:
        for choice in trans[head]:
            if len(choice) + 2 < width:
                return choice
    else:
        return head


def fline(outer, inner, guts, span=1):
    eff_width = cwidth * span + span - 1
    sguts = [f"{fheader(item, eff_width):^{eff_width}}" for item in guts]
    vertex = "|"
    return " " * pwidth + vertex + f"{' ' + outer:<{l1width}}" + vertex + f"{' ' + inner:<{l2width}}" + vertex + "|".join(sguts) + vertex


def fline_fill(outer, inner, guts, span=1):
    eff_width = cwidth * span + span - 1
    sguts = [f"{item*eff_width}" for item in guts]
    vertex = "+"
    return " " * pwidth + vertex + f"{outer*l1width}" + vertex + f"{inner*l2width}" + vertex + "+".join(sguts) + vertex


def extract_modules(winnowed):
    return sorted(list(set([entry["module"] for entry in winnowed])))


with open("stdsuite_qcng.txt", "r") as fp:
    contents = fp.readlines()

stuff = [ast.literal_eval(ln) for ln in contents]

# left margin, left outer header, left inner header, and body column widths
pwidth = 2
l1width = 25
l2width = 25
cwidth = 13  # 13 minimum for footnotes

notes = []
methods = [
"hf",
"mp2",
"mp3",
"mp4(sdq)",
"mp4",
"cisd",
"qcisd",
"qcisd(t)",
"fci",
"lccd",
"lccsd",
"ccd",
"ccsd",
"ccsd",
"ccsd+t(ccsd)",
"ccsd(t)",
"a-ccsd(t)",
"ccsdt-1a",
"ccsdt-1b",
"ccsdt-2",
"ccsdt-3",
"ccsdt",
"ccsdt(q)",
"ccsdtq",
"pbe",
"b3lyp",
"b3lyp5",
]
methods.extend([entry["method"] for entry in stuff])
methods = unique(methods)
#methods = sorted(list(set([entry["method"] for entry in stuff])))
fcaes = ["ae", "fc"]
eghs = ["energy", "gradient", "hessian"]
refs = ["rhf", "uhf", "rohf"]


#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |                                  All-Electron                                  |                                  Frozen-Core                                   |
#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |           RHF            |           UHF            |           ROHF           |           RHF            |           UHF            |           ROHF           |
#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |   E    |   G    |   H    |   E    |   G    |   H    |   E    |   G    |   H    |   E    |   G    |   H    |   E    |   G    |   H    |   E    |   G    |   H    |
#  +========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+

if False:

    lines = []
    ncol = len(fcaes) * len(refs) * len(eghs)
    lines.append(fline_fill("-", "-", ["-"] * ncol))
    lines.append(fline("", "", fcaes, span=9))
    lines.append(fline_fill(" ", " ", ["-"] * ncol))
    lines.append(fline("", "", refs*len(fcaes), span=3))
    lines.append(fline_fill(" ", " ", ["-"] * ncol))
    lines.append(fline("", "", eghs*len(fcaes)*len(refs), span=1))
    lines.append(fline_fill("=", "=", ["="] * ncol))
    
    for method in methods:
        subset1 = [entry for entry in stuff if entry["method"] == method]
        modules = extract_modules(subset1)
    
        for module in modules:
            subset2 = [entry for entry in subset1 if entry["module"] == module]
    
            line = []
    
            for corl_type in ["conv"]:
                subset3 = [entry for entry in subset2 if entry["corl_type"] == corl_type]
    
                for fcae in fcaes:
                    for ref in refs:
                        for egh in eghs:
                            subset4 = [entry for entry in subset3 if (entry["fcae"] == fcae and entry["reference"] == ref and entry["driver"] == egh)]
    
                            line.append(fcell(subset4))
            
            print_method = method if module == modules[0] else ""
            lines.append(fline(print_method, module, line))
            end_method = "-" if module == modules[-1] else " "
            #lines.append(fline_fill(end_method, "-", ["-"] * ncol))  # TOGGLE PAIR A -- for reST
    
        lines.append(fline_fill("-", "-", ["-"] * ncol))  # TOGGLE PAIR A -- for eyes
    print("\n".join(lines))


#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |                         RHF                         |                         UHF                         |                        ROHF                         |
#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |     Energy      |    Gradient     |     Hessian     |     Energy      |    Gradient     |     Hessian     |     Energy      |    Gradient     |     Hessian     |
#  +--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+
#  |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |   AE   |   FC   |
#  +========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+========+

if True:

    lines = []
    ncol = len(fcaes) * len(refs) * len(eghs)
    lines.append(fline_fill("-", "-", ["-"] * ncol))
    lines.append(fline("", "", refs, span=6))
    lines.append(fline_fill(" ", " ", ["-"] * ncol))
    lines.append(fline("", "", eghs*len(refs), span=2))
    lines.append(fline_fill(" ", " ", ["-"] * ncol))
    lines.append(fline("", "", fcaes*len(refs)*len(eghs), span=1))
    lines.append(fline_fill("=", "=", ["="] * ncol))
    
    for method in methods:
        subset1 = [entry for entry in stuff if entry["method"] == method]
        if len(subset1) == 0:
            continue

        modules = extract_modules(subset1)
        for module in modules:
            subset2 = [entry for entry in subset1 if entry["module"] == module]
    
            line = []
    
            for corl_type in ["conv"]:
                subset3 = [entry for entry in subset2 if entry["corl_type"] == corl_type]
    
                for ref in refs:
                    for egh in eghs:
                        for fcae in fcaes:
                            subset4 = [entry for entry in subset3 if (entry["fcae"] == fcae and entry["reference"] == ref and entry["driver"] == egh)]
    
                            line.append(fcell(subset4))
            
            print_method = method if module == modules[0] else ""
            lines.append(fline(print_method, module, line))
            end_method = "-" if module == modules[-1] else " "
            #lines.append(fline_fill(end_method, "-", ["-"] * ncol))  # TOGGLE PAIR A -- for reST
    
        lines.append(fline_fill("-", "-", ["-"] * ncol))  # TOGGLE PAIR A -- for eyes

    lines.append("")
    lines.append(".. rubric:: Footnotes")
    lines.append("")
    for idx, note in enumerate(notes):
        lines.append(f".. [#f{idx+10}] {note}")

    print("\n".join(lines))

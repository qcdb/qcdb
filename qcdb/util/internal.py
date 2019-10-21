import pprint
from typing import Dict

from .. import __version__

pp = pprint.PrettyPrinter(width=120)



def provenance_stamp(routine: str) -> Dict[str, str]:
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCDB's credentials for creator and version. The
    generating routine's name is passed in through `routine`.

    """
    return {'creator': 'QCDB', 'version': __version__, 'routine': routine}


def print_jobrec(label, dicary, do_print):
    """Consolidate jobrec debug printing."""

    if do_print:
        print(label + ' <<<')
        pp.pprint(dicary)
        print('>>>')

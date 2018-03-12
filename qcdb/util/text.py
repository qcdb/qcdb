

def banner(text, type=1, width=35):
    """Function to print *text* to output file in a banner of
    minimum width *width* and minimum three-line height for
    *type* = 1 or one-line height for *type* = 2. If *strNotOutfile*
    is True, function returns string rather than printing it
    to output file.

    """
    lines = text.split('\n')
    max_length = max(len(ln) for ln in lines)
    max_length = max([width, max_length])

    null = ''
    if type == 1:
        #banner = '  //' + null.center(max_length, '>') + '//\n'
        #for line in lines:
        #    banner += '  //' + line.center(max_length) + '//\n'
        #banner += '  //' + null.center(max_length, '<') + '//\n'
        banner = '  <<' + null.center(max_length, '>') + '>>\n'
        for line in lines:
            banner += '  <<' + line.center(max_length) + '>>\n'
        banner += '  <<' + null.center(max_length, '<') + '>>\n'

    if type == 2:
        banner = ''
        for line in lines:
            banner += (' ' + line + ' ').center(max_length, '=')

    return banner


def levenshtein(seq1, seq2):
    """Levenshtein distance between two strings."""

    oneago = None
    thisrow = list(range(1, len(seq2) + 1)) + [0]
    for x in range(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in range(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]


def find_approximate_string_matches(seq1, options, max_distance):
    """Compute approximate string matches from a list of options.

    Parameters
    ----------
    seq1 : str
        Target sequence to seek among `options`.
    options : list
        Valid candidates for match.
    max_distance : int
        Largest allowed difference between `seq1` and returned match.

    Returns
    -------
    list
        Items from `options` differing from `seq1` by no more than `max_distance`.

    """
    matches = []
    for seq2 in options:
        distance = levenshtein(seq1, seq2)
        if distance <= max_distance:
            matches.append(seq2)
    return matches


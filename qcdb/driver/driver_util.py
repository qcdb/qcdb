from ..util import yes, no, der0th, der1st, der2nd

def kwargs_lower(kwargs):
    """Sanitize user's `kwargs`.
    * Rebuilds and returns `kwargs` dictionary with all keys made lowercase.
    * Should be called by every function that could be called directly by the user.
    * Turns boolean-like values into actual booleans.
    * Turns values lowercase if sensible.

    """
    caseless_kwargs = {}
    for key, value in kwargs.items():
        lkey = key.lower()
        if lkey in ['subset', 'banner']:  # only kw for which case matters
            lvalue = value
        else:
            try:
                lvalue = value.lower()
            except (AttributeError, KeyError):
                lvalue = value

        if lkey in ['irrep', 'check_bsse', 'linkage', 'bsse_type']:
            caseless_kwargs[lkey] = lvalue

        elif 'dertype' in lkey:
            if der0th.match(str(lvalue)):
                caseless_kwargs[lkey] = 0
            elif der1st.match(str(lvalue)):
                caseless_kwargs[lkey] = 1
            elif der2nd.match(str(lvalue)):
                caseless_kwargs[lkey] = 2
            else:
                raise KeyError('Derivative type key ({}) was not recognized'.format(key))

        elif yes.match(str(lvalue)):
            caseless_kwargs[lkey] = True

        elif no.match(str(lvalue)):
            caseless_kwargs[lkey] = False

        else:
            caseless_kwargs[lkey] = lvalue

    return caseless_kwargs



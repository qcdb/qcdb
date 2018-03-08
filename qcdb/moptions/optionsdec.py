import sys
import uuid
from functools import wraps


def register_opts(ros):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            #accession = uuid.uuid4()
            accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
            print('<<< entering {:10} {}'.format(func.__name__, accession))
            kwargs['accession'] = accession
            ret = func(*args, **kwargs)
            print('>>> exiting  {:10} {}'.format(func.__name__, accession))
            ros.unwind_by_accession(accession)
            print(ros.print_changed())
            return ret
        return wrapper
    return decorator


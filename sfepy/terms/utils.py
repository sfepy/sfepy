import numpy as nm

def check_finiteness(data, info):
    is_finite = nm.isfinite(data)
    if not is_finite.all():
        ii = nm.where(is_finite == False)
        print ii
        print data[ii]
        msg = 'infinite %s!, see above' % info
        raise ValueError(msg)

from sfepy.base.base import *

def extend_cell_data( data, domain, rname, val = None ):
    """Extend cell data defined in a region rname to the whole domain using the
    value val, or the smallest value in data if val is None."""
    n_el = domain.shape.n_el
    if data.shape[0] == n_el: return data

    if val is None:
        if data.shape[2] > 1: # Vector.
            val = nm.amin( nm.abs( data ) )
        else: # Scalar.
            val = nm.amin( data )

    edata = nm.empty( (n_el,) + data.shape[1:], dtype = nm.float64 )
    edata.fill( val )

    region = domain.regions[rname]
    offs = region.get_cell_offsets()
    eoffs = domain.get_cell_offsets()
##     print offs
##     print eoffs
##     print domain.mat_ids_to_i_gs
##     pause()

    for group in domain.iter_groups():
        ig = group.ig
        ii = eoffs[ig]
        if ig in region.igs:
            n_cell = region.shape[ig].n_cell
            ir = offs[ig]
            edata[ii+region.cells[ig]] = data[ir:ir+n_cell]
    return edata

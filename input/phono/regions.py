def get_box_regions( dim, sizes ):
    
    if dim == 3:
        wx, wy, wz = sizes
        regions = {
            'Near' : ('nodes in (y < -%.3f)' % wy, {}),
            'Far' : ('nodes in (y > %.3f)' % wy, {}),
            'Bottom' : ('nodes in (z < -%.3f)' % wz, {}),
            'Top' : ('nodes in (z > %.3f)' % wz, {}),
            'Left' : ('nodes in (x < -%.3f)' % wx, {}),
            'Right' : ('nodes in (x > %.3f)' % wx, {}),
            'Corners' : ("""nodes in
                            ((x < -%.3f) & (y < -%.3f) & (z < -%.3f))
                          | ((x >  %.3f) & (y < -%.3f) & (z < -%.3f))
                          | ((x >  %.3f) & (y >  %.3f) & (z < -%.3f))
                          | ((x < -%.3f) & (y >  %.3f) & (z < -%.3f))
                          | ((x < -%.3f) & (y < -%.3f) & (z >  %.3f))
                          | ((x >  %.3f) & (y < -%.3f) & (z >  %.3f))
                          | ((x >  %.3f) & (y >  %.3f) & (z >  %.3f))
                          | ((x < -%.3f) & (y >  %.3f) & (z >  %.3f))
                          """ % ((wx, wy, wz) * 8), {}),
        }
    else:
        wx, wy = sizes
        regions = {
            'Bottom' : ('nodes in (y < -%.3f)' % wy, {}),
            'Top' : ('nodes in (y > %.3f)' % wy, {}),
            'Left' : ('nodes in (x < -%.3f)' % wx, {}),
            'Right' : ('nodes in (x > %.3f)' % wx, {}),
            'Corners' : ("""nodes in
                              ((x < -%.3f) & (y < -%.3f))
                            | ((x >  %.3f) & (y < -%.3f))
                            | ((x >  %.3f) & (y >  %.3f))
                            | ((x < -%.3f) & (y >  %.3f))
                            """ % ((wx, wy) * 4), {}),
        }

    return regions

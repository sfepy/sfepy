from sfepy.base.base import *
from sfepy.base.ioutils import read_dict_hdf5, write_dict_hdf5
from sfepy.homogenization.utils import iter_sym

##
# 09.07.2007, c
class Coefficients( Struct ):

    ##
    # 09.07.2007, c
    def from_file_hdf5( filename ):
        obj = Coefficients()
        obj.__dict__ = read_dict_hdf5( filename )
        for key, val in obj.__dict__.iteritems():
            if type( val ) == list:
                for ii, vv in enumerate( val ):
                    val[ii] = nm.array( vv, dtype = nm.float64 )
                    
        return obj
    from_file_hdf5 = staticmethod( from_file_hdf5 )

    ##
    # 09.07.2007, c
    def to_file_hdf5( self, filename ):
        write_dict_hdf5( filename, self.__dict__ )

    ##
    # c: 07.03.2008, r: 07.03.2008
    def _write2d( self, val ):
        fd.write( r'  \begin{equation}' )
        fd.write( '\n' )
        fd.write( r'    \left[\begin{array}{%s}' % (('c',) * val.shape[0]) )
        fd.write( '\n' )
        for ir in xrange( val.shape[1] ):
            for ic in xrange( val.shape[0] ):
                fd.write( '    ' + self.format % val[ir,ic] )
                if ic < (val.shape[0] - 1):
                    fd.write( r' & ' )
                elif ir < (val.shape[1] - 1):
                    fd.write( r' \\' )
                    fd.write( '\n' )
        fd.write( '\n' )
        fd.write( r'    \end{array}\right]' )
        fd.write( '\n' )
        fd.write( r'  \end{equation}' )
        fd.write( '\n' )
        
    ##
    # c: 07.03.2008, r: 07.03.2008
    def to_file_latex( self, filename, names, print_digits ):

        self.format = '%% %d.%df' % (print_digits + 3, print_digits)
        print self.format
        fd = open( filename, 'w' )
        fd.write( r'\begin{itemize}' )
        fd.write( '\n' )
        for key, val in self.__dict__.iteritems():
            try:
                lname = names[key]
            except:
                lname = key
            fd.write( '\item %s:' % lname )
            fd.write( '\n' )

            if val.ndim == 2:
                self._write2d( val )
        fd.write( r'\end{itemize}' )
        fd.write( '\n\n' )
        fd.close()

    ##
    # c: 07.03.2008, r: 09.04.2008
    def _save_dict( self, adict, fd, names, format ):
        for key, val in adict.iteritems():
            try:
                lname = names[key]
            except:
                lname = key
            fd.write( '%s:\n' % lname )
            if isinstance( val, dict ):
                self._save_dict( val, fd, names, format )
                fd.write( '\n' )
            elif isinstance( val, str ):
                fd.write( val + '\n' )
            elif val.ndim == 0:
                fd.write( format % val )
                fd.write( '\n' )
            elif val.ndim == 1:
                for ic in xrange( val.shape[0] ):
                    fd.write( format % val[ic] )
                    if ic < (val.shape[0] - 1):
                        fd.write( ', ' )
                    else:
                        fd.write( '\n' )
            elif val.ndim == 2:
                for ir in xrange( val.shape[0] ):
                    for ic in xrange( val.shape[1] ):
                        fd.write( format % val[ir,ic] )
                        if ic < (val.shape[1] - 1):
                            fd.write( ', ' )
                        elif ir < (val.shape[0] - 1):
                            fd.write( ';\n' )
                fd.write( '\n' )
            elif val.ndim == 3:
                for ii in xrange( val.shape[0] ):
                    fd.write( '  step %d:\n' % ii )
                    for ir in xrange( val.shape[1] ):
                        for ic in xrange( val.shape[2] ):
                            fd.write( '  ' + format % val[ii,ir,ic] )
                            if ic < (val.shape[2] - 1):
                                fd.write( ', ' )
                            elif ir < (val.shape[1] - 1):
                                fd.write( ';\n' )
                    fd.write( '\n' )
                fd.write( '\n' )
            fd.write( '\n' )

    ##
    # c: 07.03.2008, r: 07.03.2008
    def to_file_txt( self, filename, names, format ):

        fd = open( filename, 'w' )
        self._save_dict( self.__dict__, fd, names, format )
        fd.close()

    _table_vector = r"""
\begin{center}
\begin{tabular}{cc}
i & value \\
%s
\end{tabular}
\end{center}
    """

    _table_matrix_1 = r"""
\begin{center}
\begin{tabular}{cc}
ij & value \\
%s
\end{tabular}
\end{center}
    """

    _table_matrix_2 = r"""
\begin{center}
\begin{tabular}{cc}
ijkl & value \\
%s
\end{tabular}
\end{center}
    """

    _itemize = r"""
\begin{itemize}
%s
\end{itemize}
    """

    ##
    # c: 09.07.2008, r: 09.07.2008
    def _typeset( self, val, dim, style = 'table', format = '%f',
                  step = None ):
        sym = (dim + 1) * dim / 2

        mode = None
        if val.ndim == 0:
            mode = 'scalar'

        elif val.ndim == 1:
            if val.shape[0] == 1:
                mode = 'scalar'
            elif val.shape[0] == dim:
                mode = 'vector'
            elif val.shape[0] == sym:
                mode = 'matrix_t1d'
                
        elif val.ndim == 2:
            if val.shape[0] == dim:
                mode = 'matrix_2D'
            elif val.shape[0] == sym:
                mode = 'matrix_t2d'

        print mode
        out = ''
        if mode == 'scalar':
            out = format % val
        elif mode == 'vector':
            aux = ' \\\\\n'.join( [r'$_%d$ & %s' % (ir + 1, format % val[ir])
                                   for ir in xrange( dim )] )
            out = self._table_vector % aux
        elif mode == 'matrix_t1d':
            aux = ' \\\\\n'.join( [r'$_{%d%d}$ & %s' % (ir + 1, ic + 1,
                                                        format % val[ii])
                                  for ii, (ir, ic) \
                                  in enumerate( iter_sym( dim ) )] )
            out = self._table_matrix_1 % aux
        elif mode == 'matrix_2D':
            aux = ' \\\\\n'.join( [r'$_{%d%d}$ & %s' % (ir + 1, ic + 1,
                                                        format % val[ir,ic])
                                  for ir in xrange( dim )
                                  for ic in xrange( dim )] )
            out = self._table_matrix_1 % aux
        elif mode == 'matrix_t2d':
            aux = ' \\\\\n'.join( [r'$_{%d%d%d%d}$ & %s' % (irr + 1, irc + 1,
                                                            icr + 1, icc + 1,
                                                            format % val[ii,jj])
                                  for ii, (irr, irc) \
                                  in enumerate( iter_sym( dim ) )
                                  for jj, (icr, icc) \
                                  in enumerate( iter_sym( dim ) )] )
            out = self._table_matrix_2 % aux
        return out

    ##
    # c: 09.07.2008, r: 09.07.2008
    def to_latex( self, attr_name, dim, style = 'table', format = '%f',
                 step = None ):

        val = getattr( self, attr_name )
        if step is not None:
            val = val[step]

        print attr_name
        print val

        if isinstance( val, dict ):
            aux = ''
            for key, dval in val.iteritems():
                aux2 = r'\item %s : %s' % (key,
                                           self._typeset( dval, dim, style,
                                                          format, step ))
                aux = '\n'.join( (aux, aux2) )
            out = self._itemize % aux
        else:
            out = self._typeset( val, dim, style, format, step )

            
        print out
#        pause()
        return out

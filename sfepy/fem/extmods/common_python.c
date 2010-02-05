#include <stdarg.h>                                                             

#include "common.h"
#include "Python.h"

int32 g_error = 0;

static char buf[1024]; /* !!! */

#undef __FUNC__
#define __FUNC__ "output"
/*!
  @par Revision history:
  - 16.02.2004, c
*/
void output( const char *what, ... )
{
  va_list ap;

  va_start( ap, what );
  vprintf( what, ap );
  va_end( ap );
}

#undef __FUNC__
#define __FUNC__ "errput"
/*!
  @par Revision history:
  - 16.02.2004, c
  - 20.02.2004
  - 26.10.2005
*/
void errput( const char *what, ... )
{
  va_list ap;

  snprintf( buf, 1020, "**ERROR** -> %s", what );
  va_start( ap, what );
  vprintf( what, ap );
  va_end( ap );
  PyErr_SetString( PyExc_RuntimeError, "ccore error (see above)" );
  g_error++;
}

void errset( const char *msg )
{
  PyErr_SetString( PyExc_RuntimeError, msg );
  g_error++;
}

#undef __FUNC__
#define __FUNC__ "errclear"
/*!
  @par Revision history:
  - 07.06.2006, c
*/
void errclear()
{
  g_error = 0;
}

/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
static int32 al_curUsage;
static int32 al_maxUsage;
static int32 al_frags;
static AllocSpace *al_head = 0;


#undef __FUNC__
#define __FUNC__ "mem_alloc_mem"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void *mem_alloc_mem( size_t size, int lineNo, char *funName,
                    char *fileName, char *dirName )
{
  char *p;
  int32 hsize = sizeof( AllocSpaceAlign );
  int32 tsize, aux;
  float64 *endptr;
  AllocSpace *head;
  
  if (size == 0) {
    errput( "%s, %s, %s, %d: zero allocation!\n",
	    dirName, fileName, funName, lineNo );
    ERR_GotoEnd( 1 );
  }

  aux = size % sizeof( float64 );
  size += (aux) ? sizeof( float64 ) - aux : 0;
  tsize = size + hsize + sizeof( float64 );
  if ((p = (char *) PyMem_Malloc( tsize )) == 0) {
    errput( "%s, %s, %s, %d: error allocating %d bytes (current: %d).\n",
	    dirName, fileName, funName, lineNo, size, al_curUsage );
    ERR_GotoEnd( 1 );
  }
  head = (AllocSpace *) p;
  p += hsize;

  if (al_head) al_head->prev = head;
  head->next     = al_head;
  al_head         = head;
  head->prev     = 0;
  head->size     = size;
  head->id       = 1234567;
  head->lineNo   = lineNo;
  head->fileName = fileName;
  head->funName  = funName;
  head->dirName  = dirName;
  head->cookie   = AL_CookieValue;
  endptr         = (float64 *) (p + size);
  endptr[0]      = (float64) AL_CookieValue;

/*    output( "%d _> ptr: %p, head: %p, end %p\n", hsize, p, head, endptr ); */

  al_curUsage += size;
  if (al_curUsage > al_maxUsage) {
    al_maxUsage = al_curUsage;
  }
  al_frags++;

  memset( p, 0, size );

  return( (void *) p );

 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }

  return( 0 );
}

#undef __FUNC__
#define __FUNC__ "mem_free_mem"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_free_mem( void *pp, int lineNo, char *funName,
                  char *fileName, char *dirName )
{
  char *p = (char *) pp;
  int32 hsize = sizeof( AllocSpaceAlign );
  float64 *endptr;
  AllocSpace *head;
  char *phead;

  if (p == 0) return;

  phead = p - hsize;
  head = (AllocSpace *) phead;
  if (head->cookie != AL_CookieValue) {
    errput( "%s, %s, %s, %d: ptr: %p, cookie: %d\n",
	    dirName, fileName, funName, lineNo,
	    p, head->cookie );
    if (head->cookie == AL_AlreadyFreed) {
      errput( "memory was already freed!\n" );
    }
    ERR_GotoEnd( 1 );
  }
  head->cookie = AL_AlreadyFreed;

  endptr = (float64 *) (p + head->size);
  if (endptr[0] != AL_CookieValue) {
    errput( "%s %s %s %d:\n",
	    dirName, fileName, funName, lineNo );
    if (endptr[0] == AL_AlreadyFreed) {
      errput( "already freed!\n" );
    } else {
      errput( "damaged tail!\n" );
    }
    ERR_GotoEnd( 1 );
  }

  endptr[0] = (float64) AL_AlreadyFreed;

  al_curUsage -= head->size;
  al_frags--;

  if (head->prev) head->prev->next = head->next;
  else al_head = head->next;

  if (head->next) head->next->prev = head->prev;

  PyMem_Free( phead );

  return;
 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }
}

#undef __FUNC__
#define __FUNC__ "mem_checkIntegrity"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_checkIntegrity( int lineNo, char *funName,
			 char *fileName, char *dirName )
{
  char *p, *pp;
  int32 cnt, allocated;
  int32 hsize = sizeof( AllocSpaceAlign );
  float64 *endptr;
  AllocSpace *head = al_head;

  output( "checking memory integrity in\n" );
  output( "%s, %s, %s(), %d:\n",
	  dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage );
  output( "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  if (head == 0) {
    goto end_label_ok;
  }

  cnt = 0;
  allocated = 0;
  while (head) {
    p = (char *) head;

    pp = p + hsize;
    if (head->cookie != AL_CookieValue) {
      errput( "ptr: %p, ptrhead: %p, cookie: %d\n",
	      pp, p , head->cookie );
      if (head->cookie == AL_AlreadyFreed) {
	errput( "memory was already freed!\n" );
      }
      ERR_GotoEnd( 1 );
    }

    endptr = (float64 *) (pp + head->size);
    if (endptr[0] != AL_CookieValue) {
      output( "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	      head->dirName, head->fileName, head->funName, head->lineNo,
	      head->size, pp );
      if (endptr[0] == AL_AlreadyFreed) {
	errput( "already freed!\n" );
      } else {
	errput( "damaged tail!\n" );
      }
      ERR_GotoEnd( 1 );
    }
    cnt++;
    allocated += head->size;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_GotoEnd( 1 );
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_GotoEnd( 1 );
  }
  if (allocated != al_curUsage) {
    errput( "memory leak!? (%d == %d)\n", allocated, al_curUsage );
    ERR_GotoEnd( 1 );
  }

 end_label_ok:
  output( "memory OK.\n" );
  return;

 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }
}

#undef __FUNC__
#define __FUNC__ "mem_statistics"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_statistics( int lineNo, char *funName,
		     char *fileName, char *dirName )
{
  output( "%s, %s, %s(), %d: memory max: %d, current: %d\n",
	  dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage );
}

#undef __FUNC__
#define __FUNC__ "mem_print"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
  - 06.10.2005
*/
int32 mem_print( FILE *file, int32 mode )
{
  int32 cnt = 0;
  int32 hsize = sizeof( AllocSpaceAlign );
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  fprintf( file, "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  if (head == 0) {
    goto end_label_ok;
  }

  while (head) {
    p = (char *) head;
/*      fprintf( file, "%d _> head: %p, end %p\n", hsize, p, */
/*  	    p + hsize + head->size ); */
    fprintf( file, "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	    head->dirName, head->fileName, head->funName, head->lineNo,
	    head->size, p + hsize );
    cnt++;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_GotoEnd( 1 );
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_GotoEnd( 1 );
  }

 end_label_ok:
  fprintf( file, "done.\n" );

  return( RET_OK );

 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }
  return( RET_Fail );
}

#undef __FUNC__
#define __FUNC__ "mem_printSome"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
  - 06.10.2005
*/
int32 mem_printSome( FILE *file, int32 mode, int32 num )
{
  int32 cnt = 0;
  int32 hsize = sizeof( AllocSpaceAlign );
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  fprintf( file, "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  fprintf( file, "printing max: %d\n", num );
  if (head == 0) {
    goto end_label_ok;
  }

  while (head) {
    p = (char *) head;
/*      fprintf( file, "%d _> head: %p, end %p\n", hsize, p, */
/*  	    p + hsize + head->size ); */
    fprintf( file, "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	    head->dirName, head->fileName, head->funName, head->lineNo,
	    head->size, p + hsize );
    cnt++;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_GotoEnd( 1 );
    }
    if (cnt == num) break;
    head = head->next;
  }

 end_label_ok:
  fprintf( file, "done.\n" );

  return( RET_OK );

 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }
  return( RET_Fail );
}

#undef __FUNC__
#define __FUNC__ "mem_freeGarbage"
/*!
  Free all memory records.
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
int32 mem_freeGarbage()
{
  int32 cnt = 0, frags = al_frags;
  int32 hsize = sizeof( AllocSpaceAlign );
  char *p;

  output( "freeing garbage.\n" );
  while (al_head) {
    p = (char *) al_head + hsize;
    free_mem( p );
/*      output( "  %s, %s, %s, %d: size: %d, ptr: %p\n", */
/*  	    head->dirName, head->fileName, head->funName, head->lineNo, */
/*  	    head->size, p + hsize ); */
    cnt++;
/*      printf( "%d %d\n", cnt, frags ); */
    if (cnt > frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_GotoEnd( 1 );
    }
  }
  if (cnt < frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_GotoEnd( 1 );
  }

  return( RET_OK );

 end_label:
  if (ERR_Chk) {
    errput( ErrHead "error exit!\n" );
  }
  return( RET_Fail );
}

#if SFEPY_PLATFORM == 0
#include <termios.h> /* tcgetattr(), tcsetattr() */
#include <unistd.h> /* read() */
#endif

#undef __FUNC__
#define __FUNC__ "sys_getch"
/*!
  @par Revision history:
  - 21.05.2002, c
*/
int sys_getch( void )
{
  char ch = 0;
#if SFEPY_PLATFORM == 0
  if (read (STDERR_FILENO, &ch, 1) < 0) {
    return( RET_Fail );
  }
#endif
  return( ch );
}

#if SFEPY_PLATFORM == 0
static struct termios term;
#endif

#undef __FUNC__
#define __FUNC__ "sys_keyboardEnableRaw"
/*!
  Set keyboard to raw mode so getch will work

  @par Revision history:
  - 21.05.2002, c
*/
void sys_keyboardEnableRaw()
{
#if SFEPY_PLATFORM == 0
  // set to non canonical mode, echo off, ignore signals
  struct termios current;
  // save current terminal settings
  tcgetattr (STDERR_FILENO, &current);

  // set to non canonical mode, echo off, ignore signals
  term = current;
  current.c_lflag &= ~(ECHO | ICANON | IEXTEN);
  current.c_cc[VMIN] = 1;
  current.c_cc[VTIME] = 0;
  tcsetattr (STDERR_FILENO, TCSAFLUSH, &current);
#endif
}

#undef __FUNC__
#define __FUNC__ "sys_keyboardDisableRaw"
/*!
  @par Revision history:
  - 21.05.2002, c
*/
void sys_keyboardDisableRaw()
{
#if SFEPY_PLATFORM == 0
  // Restore old terminal settings
  tcsetattr (STDERR_FILENO, TCSAFLUSH, &term);
#endif
}

#undef __FUNC__
#define __FUNC__ "sys_pause()"
/*!
  @par Revision history:
  - 18.02.2005, c
  - 28.11.2005
*/
void sys_pause()
{
  sys_keyboardEnableRaw();
  if (sys_getch() == 'q') {
    sys_keyboardDisableRaw();
    exit( 1 );
  }
  sys_keyboardDisableRaw();
}

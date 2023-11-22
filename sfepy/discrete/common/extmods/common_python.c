#include <stdarg.h>

#include "common.h"

int32 g_error = 0;

#undef __FUNC__
#define __FUNC__ "output"
/*!
  @par Revision history:
  - 16.02.2004, c
*/
void output(const char *what, ...)
{
  va_list ap;

  va_start(ap, what);
  vprintf(what, ap);
  va_end(ap);
}

#undef __FUNC__
#define __FUNC__ "errput"
/*!
  @par Revision history:
  - 16.02.2004, c
  - 20.02.2004
  - 26.10.2005
*/
void errput(const char *what, ...)
{
  va_list ap;

  va_start(ap, what);
  vprintf(what, ap);
  va_end(ap);
  PyErr_SetString(PyExc_RuntimeError, "ccore error (see above)");
  g_error++;
}

void errset(const char *msg)
{
  PyErr_SetString(PyExc_RuntimeError, msg);
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
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
  }
  g_error = 0;
}

/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
static size_t al_curUsage;
static size_t al_maxUsage;
static size_t al_frags;
static AllocSpace *al_head = 0;

size_t mem_get_cur_usage(void)
{
  return al_curUsage;
}

size_t mem_get_max_usage(void)
{
  return al_maxUsage;
}

size_t mem_get_n_frags(void)
{
  return al_frags;
}

void mem_list_new(char *p, size_t size, AllocSpace *al_head,
                  int lineNo, char *funName, char *fileName, char *dirName)
{
  float64 *endptr;
  size_t hsize = sizeof(AllocSpaceAlign);
  AllocSpace *head = (AllocSpace *) (p - hsize);

  if (al_head) al_head->prev = head;
  head->next     = al_head;
  al_head        = head;
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
}

void mem_list_remove(AllocSpace *head, AllocSpace *al_head)
{
  if (head->prev) head->prev->next = head->next;
  else al_head = head->next;

  if (head->next) head->next->prev = head->prev;
}

int32 mem_check_ptr(char *p, int lineNo, char *funName,
                    char *fileName, char *dirName)
{
  int32 ret = RET_OK;
  float64 *endptr;
  size_t hsize = sizeof(AllocSpaceAlign);
  AllocSpace *head = (AllocSpace *) (p - hsize);

  if (head->cookie != AL_CookieValue) {
    errput("%s, %s, %s, %d: ptr: %p, cookie: %d\n",
           dirName, fileName, funName, lineNo,
           p, head->cookie);
    if (head->cookie == AL_AlreadyFreed) {
      errput("memory was already freed!\n");
    }
    ERR_CheckGo(ret);
  }

  endptr = (float64 *) (p + head->size);
  if (endptr[0] != AL_CookieValue) {
    errput("%s %s %s %d:\n",
           dirName, fileName, funName, lineNo);
    if (endptr[0] == AL_AlreadyFreed) {
      errput("already freed!\n");
    } else {
      errput("damaged tail!\n");
    }
    ERR_CheckGo(ret);
  }

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "mem_alloc_mem"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void *mem_alloc_mem(size_t size, int lineNo, char *funName,
                    char *fileName, char *dirName)
{
  char *p;
  size_t hsize = sizeof(AllocSpaceAlign);
  size_t tsize, aux;

  if (size == 0) {
    errput("%s, %s, %s, %d: zero allocation!\n",
           dirName, fileName, funName, lineNo);
    ERR_GotoEnd(1);
  }

  aux = size % sizeof(float64);
  size += (aux) ? sizeof(float64) - aux : 0;
  tsize = size + hsize + sizeof(float64);
  if ((p = (char *) PyMem_Malloc(tsize)) == 0) {
    errput("%s, %s, %s, %d: error allocating %zu bytes (current: %zu).\n",
           dirName, fileName, funName, lineNo, size, al_curUsage);
    ERR_GotoEnd(1);
  }
  p += hsize;

  mem_list_new(p, size, al_head, lineNo, funName, fileName, dirName);

  al_curUsage += size;
  if (al_curUsage > al_maxUsage) {
    al_maxUsage = al_curUsage;
  }
  al_frags++;

  memset(p, 0, size);

  return((void *) p);

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }

  return(0);
}

#undef __FUNC__
#define __FUNC__ "mem_realloc_mem"
void *mem_realloc_mem(void *pp, size_t size, int lineNo, char *funName,
                      char *fileName, char *dirName)
{
  char *p = (char *) pp;
  size_t hsize = sizeof(AllocSpaceAlign);
  size_t tsize, aux;
  float64 *endptr;
  AllocSpace *head;
  char *phead;

  if (p == 0) return(0);

  if (size == 0) {
    errput("%s, %s, %s, %d: zero allocation!\n",
           dirName, fileName, funName, lineNo);
    ERR_GotoEnd(1);
  }

  // 1. almost as mem_free_mem().
  mem_check_ptr(p, lineNo, funName, fileName, dirName);
  if (ERR_Chk) {
    ERR_GotoEnd(1);
  }

  phead = p - hsize;
  head = (AllocSpace *) phead;
  head->cookie = AL_AlreadyFreed;

  endptr = (float64 *) (p + head->size);
  endptr[0] = (float64) AL_AlreadyFreed;

  al_curUsage -= head->size;
  al_frags--;
  mem_list_remove(head, al_head);

  // 2. realloc.
  aux = size % sizeof(float64);
  size += (aux) ? sizeof(float64) - aux : 0;
  tsize = size + hsize + sizeof(float64);
  if ((p = (char *) PyMem_Realloc(phead, tsize)) == 0) {
    errput("%s, %s, %s, %d: error re-allocating to %zu bytes (current: %zu).\n",
           dirName, fileName, funName, lineNo, size, al_curUsage);
    ERR_GotoEnd(1);
  }

  // 3. almost as mem_alloc_mem().
  p += hsize;
  mem_list_new(p, size, al_head, lineNo, funName, fileName, dirName);

  al_curUsage += size;
  if (al_curUsage > al_maxUsage) {
    al_maxUsage = al_curUsage;
  }
  al_frags++;

  return((void *) p);

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }

  return(0);
}

#undef __FUNC__
#define __FUNC__ "mem_free_mem"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_free_mem(void *pp, int lineNo, char *funName,
                  char *fileName, char *dirName)
{
  char *p = (char *) pp;
  size_t hsize = sizeof(AllocSpaceAlign);
  float64 *endptr;
  AllocSpace *head;
  char *phead;

  if (p == 0) return;

  mem_check_ptr(p, lineNo, funName, fileName, dirName);
  if (ERR_Chk) {
    ERR_GotoEnd(1);
  }

  phead = p - hsize;
  head = (AllocSpace *) phead;
  head->cookie = AL_AlreadyFreed;

  endptr = (float64 *) (p + head->size);
  endptr[0] = (float64) AL_AlreadyFreed;

  al_curUsage -= head->size;
  al_frags--;

  mem_list_remove(head, al_head);

  PyMem_Free(phead);

  return;

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }
}

#undef __FUNC__
#define __FUNC__ "pyalloc"
void *pyalloc(size_t size)
{
  return mem_alloc_mem(size, __LINE__, __FUNC__, __FILE__, __SDIR__);
}

#undef __FUNC__
#define __FUNC__ "pyfree"
void pyfree(void *pp)
{
  return mem_free_mem(pp, __LINE__, __FUNC__, __FILE__, __SDIR__);
}

#undef __FUNC__
#define __FUNC__ "mem_checkIntegrity"
/*!
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_checkIntegrity(int lineNo, char *funName,
                        char *fileName, char *dirName)
{
  char *p, *pp;
  size_t cnt, allocated;
  size_t hsize = sizeof(AllocSpaceAlign);
  float64 *endptr;
  AllocSpace *head = al_head;

  output("checking memory integrity in\n");
  output("%s, %s, %s(), %d:\n",
         dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage);
  output("allocated memory: %zu records, usage: %zu, max: %zu\n",
         al_frags, al_curUsage, al_maxUsage);
  if (head == 0) {
    goto end_label_ok;
  }

  cnt = 0;
  allocated = 0;
  while (head) {
    p = (char *) head;

    pp = p + hsize;
    if (head->cookie != AL_CookieValue) {
      errput("ptr: %p, ptrhead: %p, cookie: %d\n",
             pp, p , head->cookie);
      if (head->cookie == AL_AlreadyFreed) {
        errput("memory was already freed!\n");
      }
      ERR_GotoEnd(1);
    }

    endptr = (float64 *) (pp + head->size);
    if (endptr[0] != AL_CookieValue) {
      output("  %s, %s, %s, %d: size: %zu, ptr: %p\n",
             head->dirName, head->fileName, head->funName, head->lineNo,
             head->size, pp);
      if (endptr[0] == AL_AlreadyFreed) {
        errput("already freed!\n");
      } else {
        errput("damaged tail!\n");
      }
      ERR_GotoEnd(1);
    }
    cnt++;
    allocated += head->size;
    if (cnt > al_frags) {
      errput("damaged allocation record (overrun)!\n");
      ERR_GotoEnd(1);
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput("damaged allocation record (underrun)!\n");
    ERR_GotoEnd(1);
  }
  if (allocated != al_curUsage) {
    errput("memory leak!? (%zu == %zu)\n", allocated, al_curUsage);
    ERR_GotoEnd(1);
  }

 end_label_ok:
  output("memory OK.\n");
  return;

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }
}

#undef __FUNC__
#define __FUNC__ "mem_statistics"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
*/
void mem_statistics(int lineNo, char *funName,
                    char *fileName, char *dirName)
{
  output("%s, %s, %s(), %d: memory max: %zu, current: %zu\n",
         dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage);
}

#undef __FUNC__
#define __FUNC__ "mem_print"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
  - 06.10.2005
*/
int32 mem_print(FILE *file, int32 mode)
{
  size_t cnt = 0;
  size_t hsize = sizeof(AllocSpaceAlign);
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  fprintf(file, "allocated memory: %zu records, usage: %zu, max: %zu\n",
          al_frags, al_curUsage, al_maxUsage);
  if (head == 0) {
    goto end_label_ok;
  }

  while (head) {
    p = (char *) head;
    fprintf(file, "  %s, %s, %s, %d: size: %zu, ptr: %p\n",
            head->dirName, head->fileName, head->funName, head->lineNo,
            head->size, p + hsize);
    cnt++;
    if (cnt > al_frags) {
      errput("damaged allocation record (overrun)!\n");
      ERR_GotoEnd(1);
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput("damaged allocation record (underrun)!\n");
    ERR_GotoEnd(1);
  }

 end_label_ok:
  fprintf(file, "done.\n");

  return(RET_OK);

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }
  return(RET_Fail);
}

#undef __FUNC__
#define __FUNC__ "mem_printSome"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 17.02.2005, from rcfem2
  - 06.10.2005
*/
int32 mem_printSome(FILE *file, int32 mode, int32 num)
{
  size_t cnt = 0;
  size_t hsize = sizeof(AllocSpaceAlign);
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  fprintf(file,
          "allocated memory: %zu records, usage: %zu, max: %zu\n",
          al_frags, al_curUsage, al_maxUsage);
  fprintf(file, "printing max: %d\n", num);
  if (head == 0) {
    goto end_label_ok;
  }

  while (head) {
    p = (char *) head;
    fprintf(file, "  %s, %s, %s, %d: size: %zu, ptr: %p\n",
            head->dirName, head->fileName, head->funName, head->lineNo,
            head->size, p + hsize);
    cnt++;
    if (cnt > al_frags) {
      errput("damaged allocation record (overrun)!\n");
      ERR_GotoEnd(1);
    }
    if (cnt == (size_t)num) break;
    head = head->next;
  }

 end_label_ok:
  fprintf(file, "done.\n");

  return(RET_OK);

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }
  return(RET_Fail);
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
  size_t cnt = 0, frags = al_frags;
  size_t hsize = sizeof(AllocSpaceAlign);
  char *p;

  output("freeing garbage.\n");
  while (al_head) {
    p = (char *) al_head + hsize;
    free_mem(p);
    cnt++;
    if (cnt > frags) {
      errput("damaged allocation record (overrun)!\n");
      ERR_GotoEnd(1);
    }
  }
  if (cnt < frags) {
    errput("damaged allocation record (underrun)!\n");
    ERR_GotoEnd(1);
  }

  return(RET_OK);

 end_label:
  if (ERR_Chk) {
    errput(ErrHead "error exit!\n");
  }
  return(RET_Fail);
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
int sys_getch(void)
{
  char ch = 0;
#if SFEPY_PLATFORM == 0
  if (read (STDERR_FILENO, &ch, 1) < 0) {
    return(RET_Fail);
  }
#endif
  return(ch);
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
    exit(1);
  }
  sys_keyboardDisableRaw();
}

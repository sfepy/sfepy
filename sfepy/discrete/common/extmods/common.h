#ifndef _COMMON_H_
#define _COMMON_H_ 1

#ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

#include "types.h"
#include "version.h"

typedef enum ReturnStatus {
  RET_OK,
  RET_Fail
} ReturnStatus;

void output(const char *what, ...);
void errput(const char *what, ...);
void errset(const char *msg);
void errclear(void);

#define AL_CookieValue   0xf0e0d0c9
#define AL_AlreadyFreed  0x0f0e0d9c

/*!
  @par Revision history:
  - 25.05.2001, c
*/
typedef struct _AllocSpace {
    unsigned long   size;
    int             id;
    int             lineNo;
    char            *fileName;
    char            *funName;
    char            *dirName;
    unsigned long   cookie;
    struct _AllocSpace *next,*prev;
} AllocSpace;

#define AL_HeaderDoubles      sizeof(AllocSpace)/sizeof(double)+1

/*!
  This union is used to insure that the block passed to the user is
  aligned on a double boundary

  @par Revision history:
  - 25.05.2001, c
*/
typedef union {
    AllocSpace sp;
    double  v[AL_HeaderDoubles];
} AllocSpaceAlign;

/*! @enum AllocMode
  @par Revision history:
  - 15.04.2001, c
*/
typedef enum AllocMode {
  AL_Alloc, AL_Free, AL_Realloc
} AllocMode;

size_t mem_get_cur_usage(void);
size_t mem_get_max_usage(void);
size_t mem_get_n_frags(void);
void mem_list_new(char *p, size_t size, AllocSpace *al_head,
                  int lineNo, char *funName, char *fileName, char *dirName);
void mem_list_remove(AllocSpace *head, AllocSpace *al_head);
int32 mem_check_ptr(char *p, int lineNo, char *funName,
                    char *fileName, char *dirName);
void *mem_alloc_mem(size_t size, int lineNo, char *funName,
                    char *fileName, char *dirName);
void *mem_realloc_mem(void *pp, size_t size, int lineNo, char *funName,
                      char *fileName, char *dirName);
void mem_free_mem(void *pp, int lineNo, char *funName,
                  char *fileName, char *dirName);
void *pyalloc(size_t size);
void pyfree(void *pp);
void mem_checkIntegrity(int lineNo, char *funName,
                        char *fileName, char *dirName);
void mem_statistics(int lineNo, char *funName,
                    char *fileName, char *dirName);
int32 mem_print(FILE *file, int32 mode);
int32 mem_printSome(FILE *file, int32 mode, int32 num);
int32 mem_freeGarbage(void);

int sys_getch(void);
void sys_keyboardEnableRaw(void);
void sys_keyboardDisableRaw(void);
void sys_pause(void);

#ifndef __SDIR__
  #define __SDIR__ ""
#endif

/*!
  @par Revision history:
  - 06.03.2003, c
*/
#define alloc_mem(Type, num) \
  (Type *) mem_alloc_mem((num) * sizeof(Type), \
			 __LINE__, __FUNC__, __FILE__, __SDIR__)
#define realloc_mem(p, Type, num)\
  (Type *) mem_realloc_mem(p, (num) * sizeof(Type), \
                           __LINE__, __FUNC__, __FILE__, __SDIR__)
#define free_mem(p) do {\
    mem_free_mem(p, __LINE__, __FUNC__, __FILE__, __SDIR__); } while (0)

#define print_mem_stats() \
  (mem_statistics(__LINE__, __FUNC__, __FILE__, __SDIR__))

#define check_memory_integrity() \
  (mem_checkIntegrity(__LINE__, __FUNC__, __FILE__, __SDIR__))

/*!
  @par Revision history:
  - 11.03.2003, c
  - 28.03.2003
  - 08.04.2003
  - 26.10.2005
*/
#define ERR_CheckGo(ret) do {\
  if (g_error != 0) {\
    (ret) = RET_Fail;\
    goto end_label;\
  }\
} while (0)
#define ERR_GotoEnd(i) do { g_error = (i); goto end_label; } while (0)
#define ERR_Chk (g_error != 0)
#define ERR_Clear (g_error = 0)
#define ErrHead __FUNC__ "(): "
extern int32 g_error;

#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Min(a,b) (((a) < (b)) ? (a) : (b))

#endif /* !SIC_COMMON_H */

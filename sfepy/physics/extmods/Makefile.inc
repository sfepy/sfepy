FEMLOCDIR := sfepy/fem/extmods
LOCDIR := sfepy/physics/extmods

LINKFLAGS := $(shell $(RUNCONFIG) link_flags)

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

SRC_SRCSWIG := $(wildcard $(LOCDIR)/*.i)
SRC_OBJSWIGC := $(SRC_SRCSWIG:.i=_wrap.c)
SRC_OBJSWIG := $(SRC_OBJSWIGC:.c=.o)
SRC_OBJSWIGPY := $(SRC_SRCSWIG:.i=.py)

DEPSF := dft.c dft_wrap.c
DEPSF := $(addprefix $(LOCDIR)/,$(DEPSF:.c=.o))
$(LOCDIR)/dft_wrap.c:
$(LOCDIR)/_dft.so: $(DEPSF)
	@echo $(DEPSF)
	$(CC) $(LINKFLAGS) $(DEPSF) $(SWIG_LIB) -o $@

SRC_LIBSWIG += $(LOCDIR)/_dft.so

SRC_SRCC := $(wildcard $(LOCDIR)/*.c)
SRC_SRCC := $(filter-out $(SRC_OBJSWIGC), $(SRC_SRCC))
SRC_OBJC := $(SRC_SRCC:.c=.o)
SRC_HDRC := $(wildcard $(LOCDIR)/*.h)

#######

SRCCFILES += $(SRC_OBJSWIGC) $(SRC_SRCC)
HDRCFILES +=
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(SRC_OBJSWIGC) $(SRC_OBJSWIG) $(SRC_OBJSWIGPY) $(SRC_LIBSWIG) $(SRC_OBJC)


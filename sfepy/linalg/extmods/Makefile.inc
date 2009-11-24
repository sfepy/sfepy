FEMLOCDIR := sfepy/fem/extmods
LOCDIR := sfepy/linalg/extmods

LINKFLAGS := $(shell $(RUNCONFIG) link_flags)

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

SRC_SRCSWIG := $(wildcard $(LOCDIR)/*.i)
SRC_OBJSWIGC := $(SRC_SRCSWIG:.i=_wrap.c)
SRC_OBJSWIG := $(SRC_OBJSWIGC:.c=.o)
SRC_OBJSWIGPY := $(SRC_SRCSWIG:.i=.py)

FEMDEPSG = common_python.c
DEPSG := rcm.c rcm_wrap.c
DEPSG := $(addprefix $(FEMLOCDIR)/,$(FEMDEPSG:.c=.o)) $(addprefix $(LOCDIR)/,$(DEPSG:.c=.o))
$(LOCDIR)/rcm_wrap.c : $(FEMLOCDIR)/common.i $(FEMLOCDIR)/array.i
$(LOCDIR)/_rcm.so: $(DEPSG)
	@echo $(DEPSG)
	$(CC) $(LINKFLAGS) $(DEPSG) $(SWIG_LIB) -o $@

SRC_LIBSWIG += $(LOCDIR)/_rcm.so

SRC_SRCC := $(wildcard $(LOCDIR)/*.c)
SRC_SRCC := $(filter-out $(SRC_OBJSWIGC), $(SRC_SRCC))
SRC_OBJC := $(SRC_SRCC:.c=.o)
SRC_HDRC := $(wildcard $(LOCDIR)/*.h)


#######

SRCCFILES += $(SRC_OBJSWIGC) $(SRC_SRCC)
HDRCFILES += 
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(SRC_OBJSWIGC) $(SRC_OBJSWIG) $(SRC_OBJSWIGPY) $(SRC_LIBSWIG) $(SRC_OBJC)

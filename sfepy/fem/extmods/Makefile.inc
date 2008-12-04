LOCDIR := sfepy/fem/extmods

LINKFLAGS := $(shell script/config.py link_flags)

SRC_SRCPY := $(shell ls $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

SRC_SRCSWIG := $(shell ls $(LOCDIR)/*.i)
#SRC_SRCSWIG := $(LOCDIR)/meshutils.i
SRC_OBJSWIGC := $(SRC_SRCSWIG:.i=_wrap.c)
SRC_OBJSWIG := $(SRC_OBJSWIGC:.c=.o)
SRC_OBJSWIGPY := $(SRC_SRCSWIG:.i=.py)

DEPSB := geomtrans.c meshutils.c meshutils_wrap.c common_python.c sort.c
DEPSB := $(addprefix $(LOCDIR)/,$(DEPSB:.c=.o))
$(LOCDIR)/meshutils_wrap.c : $(LOCDIR)/common.i $(LOCDIR)/array.i
$(LOCDIR)/_meshutils.so: $(DEPSB)
	@echo $(DEPSB)
	$(CC) $(LINKFLAGS) $(DEPSB) $(SWIG_LIB) -o $@

DEPSC := fmfield.c geometry.c geometry_wrap.c geommech.c common_python.c
DEPSC := $(addprefix $(LOCDIR)/,$(DEPSC:.c=.o))
$(LOCDIR)/geometry_wrap.c : $(LOCDIR)/common.i $(LOCDIR)/array.i $(LOCDIR)/fmfield.i
$(LOCDIR)/_geometry.so: $(DEPSC)
	@echo $(DEPSC)
	$(CC) $(LINKFLAGS) $(DEPSC) $(SWIG_LIB) -o $@

DEPSE := fmfield.c fem.c fem_wrap.c sort.c common_python.c
DEPSE := $(addprefix $(LOCDIR)/,$(DEPSE:.c=.o))
$(LOCDIR)/fem_wrap.c : $(LOCDIR)/common.i $(LOCDIR)/array.i $(LOCDIR)/fmfield.i
$(LOCDIR)/_fem.so: $(DEPSE)
	@echo $(DEPSE)
	$(CC) $(LINKFLAGS) $(DEPSE) $(SWIG_LIB) -o $@

SRC_LIBSWIG += $(LOCDIR)/_meshutils.so $(LOCDIR)/_geometry.so $(LOCDIR)/_fem.so

SRC_SRCC := $(shell ls $(LOCDIR)/*.c)
SRC_SRCC := $(filter-out $(SRC_OBJSWIGC), $(SRC_SRCC))
SRC_OBJC := $(SRC_SRCC:.c=.o)
SRC_HDRC := $(shell ls $(LOCDIR)/*.h)


#######

SRCCFILES += $(SRC_OBJSWIGC) $(SRC_SRCC)
HDRCFILES += 
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(SRC_OBJSWIGC) $(SRC_OBJSWIG) $(SRC_OBJSWIGPY) $(SRC_LIBSWIG) $(SRC_OBJC)

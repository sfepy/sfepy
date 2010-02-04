FEMLOCDIR := sfepy/fem/extmods
LOCDIR := sfepy/terms/extmods

LINKFLAGS := $(shell $(RUNCONFIG) link_flags)

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

SRC_SRCSWIG := $(wildcard $(LOCDIR)/*.i)
SRC_OBJSWIGC := $(SRC_SRCSWIG:.i=_wrap.c)
SRC_OBJSWIG := $(SRC_OBJSWIGC:.c=.o)
SRC_OBJSWIGPY := $(SRC_SRCSWIG:.i=.py)

FEMDEPSD = common_python.c fmfield.c geommech.c
DEPSD := formSDCC.c terms.c termsBasic.c termsAdjointNavierStokes.c termsElectric.c termsHDPM.c termsMass.c termsNavierStokes.c termsBiot.c termsLaplace.c termsLinElasticity.c termsHyperElasticity.c termsPiezo.c termsSurface.c termsVolume.c terms_wrap.c termsAcoustic.c
ifndef ISOPT
  DEPSD := $(filter-out termsAdjointNavierStokes.c, $(DEPSD))
endif
DEPSD := $(addprefix $(FEMLOCDIR)/,$(FEMDEPSD:.c=.o)) $(addprefix $(LOCDIR)/,$(DEPSD:.c=.o))
$(LOCDIR)/terms_wrap.c : $(FEMLOCDIR)/common.i $(FEMLOCDIR)/array.i $(FEMLOCDIR)/fmfield.i
$(LOCDIR)/_terms.so: $(DEPSD)
	@echo $(DEPSD)
	$(CC) $(LINKFLAGS) $(DEPSD) $(SWIG_LIB) -o $@


SRC_LIBSWIG += $(LOCDIR)/_terms.so

SRC_SRCC := $(wildcard $(LOCDIR)/*.c)
SRC_SRCC := $(filter-out $(SRC_OBJSWIGC), $(SRC_SRCC))
SRC_OBJC := $(SRC_SRCC:.c=.o)
SRC_HDRC := $(wildcard $(LOCDIR)/*.h)

#######

SRCCFILES += $(SRC_OBJSWIGC) $(SRC_SRCC)
HDRCFILES += 
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(SRC_OBJSWIGC) $(SRC_OBJSWIG) $(SRC_OBJSWIGPY) $(SRC_LIBSWIG) $(SRC_OBJC)

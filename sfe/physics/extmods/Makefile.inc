# c: 21.03.2008, r: 21.03.2008
LOCDIR := sfe/physics/extmods

SRC_SRCPY := $(shell ls $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

DEPSA := rdirac.f
DEPSA := $(addprefix $(LOCDIR)/,$(DEPSA))
$(LOCDIR)/rdirac.so: $(DEPSA)
	@echo $(DEPSA)
	f2py --fcompiler=gnu95 --f77flags="-W -Wall -fdefault-real-8" -c -m rdirac $(DEPSA)
	mv rdirac.so $@

SRC_LIBSWIG += $(LOCDIR)/rdirac.so

#######

SRCCFILES +=
HDRCFILES += 
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(LOCDIR)/rdirac.so

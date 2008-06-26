# c: 21.03.2008, r: 21.03.2008
LOCDIR := sfepy/physics/extmods

SRC_SRCPY := $(shell ls $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

DEPSA := dft.f
DEPSA := $(addprefix $(LOCDIR)/,$(DEPSA))
$(LOCDIR)/dft.so: $(DEPSA)
	@echo $(DEPSA)
	f2py --fcompiler=gnu95 --f77flags="-W -Wall -fdefault-real-8" -c -m dft $(DEPSA)
	mv dft.so $@

SRC_LIBSWIG += $(LOCDIR)/dft.so

#######

SRCCFILES +=
HDRCFILES +=
SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY) $(LOCDIR)/dft.so

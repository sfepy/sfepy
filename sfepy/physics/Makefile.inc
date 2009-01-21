# c: 21.03.2008, r: 21.03.2008
LOCDIR := sfepy/physics

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

#######

SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY)

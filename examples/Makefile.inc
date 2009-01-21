# c: 05.05.2008
LOCDIR := examples

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

#######

SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY)

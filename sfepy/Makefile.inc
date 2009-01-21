# 01.08.2006, c
# last revision: 11.01.2007
LOCDIR := sfepy

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

#######

SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY)

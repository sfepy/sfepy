# 19.09.2006, c
# last revision: 14.11.2007 
LOCDIR := tests

SRC_SRCPY := $(wildcard $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

#######

SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY)

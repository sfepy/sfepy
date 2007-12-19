# 09.03.2007, c
# last revision: 09.03.2007
LOCDIR := sfe/homogenization

SRC_SRCPY := $(shell ls $(LOCDIR)/*.py)
SRC_OBJPY := $(SRC_SRCPY:.py=.pyc)

#######

SRCPYFILES += $(SRC_SRCPY)
CLEANFILES += $(SRC_OBJPY)

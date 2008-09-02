# 14.12.2004, c
# last revision: 12.08.2008
VERSION := 00.50.00
PROJECTNAME := sfepy

############### Edit here. #######################################

C++         := gcc
CC          := gcc
SWIG        := swig
DATE        := date +%Y_%m_%d

CARCHFLAGS   := -Wall -c
CARCHOUT     := -o

#DEBUG_FLAGS := -DDEBUG_FMF -DDEBUG_MESH
DEBUG_FLAGS := -DDEBUG_FMF
#DEBUG_FLAGS :=

################ Do not edit below! ##############################

#  OPTFLAGS     := -g -pg -fPIC -DPIC
OPTFLAGS     := $(shell script/config.py opt_flags)
PYVER := $(shell script/config.py python_version)
ARCHLIB := $(shell script/config.py archlib)
NUMPYINCLUDE := $(shell script/config.py numpy_include)

PYTHON_INCL  := -I/usr/include/python$(PYVER) -I$(NUMPYINCLUDE)
#  SWIG_LIB     := -lswigpy

EXT_INCL     := $(PYTHON_INCL)

###############

ISRELEASE := 1
MODULES := eldesc examples input sfepy sfepy/base sfepy/fem sfepy/fem/extmods sfepy/homogenization sfepy/solvers sfepy/terms sfepy/terms/extmods sfepy/physics sfepy/physics/extmods tests
ifndef ISRELEASE
  MODULES += sfepy/optimize
else
  DEBUG_FLAGS :=
endif
VERSIONH := sfepy/fem/extmods/version.h
ALLTARGETS := version modules

CUR_DIR := $(shell pwd)

DISTFILES_TOP := btrace_python Makefile DIARY VERSION findSurf.py init_sfepy.py shaper.py test.mesh gen genhtml genDocs.py genPerMesh.py homogen.py extractor.py plotPerfusionCoefs.py runTests.py simple.py schroedinger.py eigen.py site_cfg_template.py TODO INSTALL README
RELDISTFILES_TOP := btrace_python Makefile VERSION init_sfepy.py extractor.py findSurf.py gen genhtml genDocs.py genPerMesh.py runTests.py simple.py schroedinger.py eigen.py site_cfg_template.py convert.py INSTALL README RELEASE_NOTES.txt
SUBDIRS = database doc eldesc examples input script sfepy tests
RELSUBDIRS = database doc eldesc examples input script sfepy tests geom
DATADIRS := database
DATADISTDIR := $(PROJECTNAME)-data-$(shell $(DATE))
DISTDIR := $(PROJECTNAME)-$(VERSION)
RELDISTDIR := $(PROJECTNAME)-release-$(VERSION)
BKUPDIR = 00backup-$(notdir $(PROJECTNAME))

####### Default rule

all:	$(ALLTARGETS)
	echo $(DATADIR)
	@echo All done.

##################################################################

INCL         :=
SRCPYFILES   := 
SRCCFILES    :=
HDRCFILES    :=
CLEANFILES   :=
SRC_LIBSWIG  :=

# include the description for each module
include $(patsubst %,%/Makefile.inc,$(MODULES))
INCL += $(patsubst %,-I%,$(MODULES))

CFLAGS := $(OPTFLAGS) -D__SDIR__='"${LOCDIR}"' ${DEBUG_FLAGS} ${INCL} ${EXT_INCL} $(CARCHFLAGS)

ifdef ISRELEASE
  CFLAGS += -DISRELEASE
  SWIGFLAGS := -DISRELEASE
endif

####### Implicit rules

%_wrap.c: %.i
#	$(SWIG) -noruntime -python $(INCL)  ${EXT_INCL} -o $@ $<
	$(SWIG) -python $(INCL)  ${EXT_INCL} $(SWIGFLAGS) -o $@ $<

%_wrap.o: %_wrap.c
	$(CC) $(CFLAGS) $< $(CARCHOUT) $@

%.o : %.c
	$(CC) $(CFLAGS) $< $(CARCHOUT) $@

####### Build rules

.PHONY : tags version dist reldist htmldocs save backup clean

modules: sfepy/fem/extmods/version.h $(SRC_LIBSWIG)
	@echo Python modules done.
	@echo ""

#$(SRC_LIBSWIG): $(SRC_OBJSWIG) $(SRC_OBJC)
#	@echo $(SRC_SRCSWIG)
#	@echo $(SRC_OBJSWIGC)
#	@echo $(SRC_OBJSWIG)
#	@echo $(SRC_OBJSWIGPY)
#	@echo $(SRC_LIBSWIG)
#	@echo $(SRC_OBJC)
#	$(CC) -shared -fPIC -DPIC $< $(SRC_OBJC) $(SWIG_LIB) -o $@
#
$(VERSIONH) : Makefile
	sed "s|^\(#define VERSION\) \".*\"|\1 \"$(shell cat VERSION)\"|;" $(VERSIONH).in > $(VERSIONH)

clean:
	-rm -f *.o *.bak *~ *% *tgz #*
	-rm -f $(CLEANFILES)

tags: clear_tags c_tags python_tags

clear_tags:
	-rm -f TAGS

c_tags:
	-etags -a $(SRCCFILES) $(HDRCFILES)

python_tags:
	-etags -a $(SRCPYFILES)

version:
	echo $(VERSION)-$(shell hg identify -i) > 'VERSION'

dist: version
	-mkdir $(DISTDIR)
	rm -rf $(DISTDIR)/*
	for i in $(DISTFILES_TOP); do cp -fpd $$i $(DISTDIR)/$$i; done
	@for i in $(SUBDIRS); do \
          $(MAKE) -C $$i -f Makefile.dist dist DISTDIR=${CUR_DIR}/$(DISTDIR); \
	done
	tar cf $(DISTDIR).tar $(DISTDIR)
	gzip -f --best $(DISTDIR).tar
	mv $(DISTDIR).tar.gz $(DISTDIR).tgz
	rm -rf $(DISTDIR)

reldist: version
	-./gen
	-mkdir $(RELDISTDIR)
	rm -rf $(RELDISTDIR)/*
	for i in $(RELDISTFILES_TOP); do cp -fpd $$i $(RELDISTDIR)/$$i; done
	@for i in $(RELSUBDIRS); do \
          $(MAKE) -C $$i -f Makefile.dist reldist DISTDIR=${CUR_DIR}/$(RELDISTDIR); \
	done
	mkdir $(RELDISTDIR)/output-tests
	mkdir $(RELDISTDIR)/tmp
	sed "s|ISRELEASE \:\=|ISRELEASE \:\= 1|;" $(RELDISTDIR)/Makefile > $(RELDISTDIR)/Makefile2
	mv -f $(RELDISTDIR)/Makefile2 $(RELDISTDIR)/Makefile
	tar cf $(RELDISTDIR).tar $(RELDISTDIR)
	gzip -f --best $(RELDISTDIR).tar
	mv -f $(RELDISTDIR).tar.gz $(RELDISTDIR).tgz
	rm -rf $(RELDISTDIR)

databackup:
	-mkdir $(DATADISTDIR)
	rm -rf $(DATADISTDIR)/*
	for i in $(DATADIRS); do cp -fa $$i $(DATADISTDIR)/$$i; done
	tar cf $(DATADISTDIR).tar $(DATADISTDIR)
	gzip -f --best $(DATADISTDIR).tar
	mv $(DATADISTDIR).tar.gz $(DATADISTDIR).tgz
	rm -rf $(DATADISTDIR)
	mv -f $(DATADISTDIR).tgz $(BKUPDIR)/$(DATADISTDIR).tgz

backup: dist
	-mkdir $(BKUPDIR)
	mv -f $(DISTDIR).tgz $(BKUPDIR)/$(DISTDIR).tgz

save: backup
	-mount /mnt/floppy
	cp -f  $(BKUPDIR)/$(DISTDIR).tgz /mnt/floppy/
	umount /mnt/floppy

htmldocs:
	-rm -rf doc/html
	-rm -f doc/doxygenrc
	-mkdir doc/aux
	-pythfilter.py . doc/aux/
	sed "s|^\(PROJECT_NUMBER         = \)NUMBER|\1$(VERSION)|;"\
	doc/doxygen.config > doc/doxygenrc
	doxygen doc/doxygenrc
	-rm -rf doc/aux

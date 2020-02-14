.DEFAULT_GOAL := all_nodebug
#include makefile_paths
# https://www.gnu.org/software/make/manual/html_node/File-Name-Functions.html  , nie przydały się...
GLOBAL_INSTALL_PARENT_DIR = /tmp/15-YADE/23-yade-install-bin/
GLOBAL_BUILD_PARENT_DIR   = /tmp/15-YADE/24-yade-kompilacja/
#GLOBAL_INSTALL_PARENT_DIR = /home/janek/15-YADE/90-Magazyn/06-intall/
#GLOBAL_BUILD_PARENT_DIR   = /home/janek/15-YADE/90-Magazyn/05-build/

SOURCE_DIR                := $(shell pwd)
# yade-gitlab-fork , yade-gitlab-fork-examples , yade-qm-trunk
SOURCE_PARENT_DIR         := $(shell basename $(SOURCE_DIR))
# useful only if there are more dirs with yade development
# in my case I have here: 00-yade-truk , 01-yade+FD-CFD , 02-yade+QM , 03-yade-molecular-dynamics
SOURCE_UPPARENT_DIR       := $(shell basename $(shell dirname $(SOURCE_DIR)))

BUILD_PARENT_DIR          := $(GLOBAL_BUILD_PARENT_DIR)$(SOURCE_UPPARENT_DIR)
BUILD_WORK_DIR            := build-$(SOURCE_PARENT_DIR)
INSTALL_DIR               := $(GLOBAL_INSTALL_PARENT_DIR)$(SOURCE_UPPARENT_DIR)/$(SOURCE_PARENT_DIR)


VERSION_GIT=$(shell git log -n1 --pretty=oneline | cut -c1-7)
VERSION_DATE=$(shell git log -n1 --pretty=fuller --date=iso | grep AuthorDate | cut -c13-22)
THE_VERSION=""
THE_VERSION=$(shell grep -s revision $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR)/config.py |  cut -c11-32 | sed -e "s/'//g" )
ifeq ($(THE_VERSION),)
THE_VERSION=$(VERSION_DATE).git-$(VERSION_GIT)
endif

# `make all` target by default shows a message when compilation is finished. I coudn't get it to work for other targets.
# `make all MESG=OFF` will disable this.
# the message program (below) is shell invocation of `mesr`. You can replace it with `xmessage` or anything you want.
MESG:=ON
GOAL_HERE::=$(MAKECMDGOALS)
# make all MESG=ON will print a message at the end of compilation (here I use a message in polish)
MAYBE_MESSAGE=$(if $(subst $@,,$(GOAL_HERE)),,$(if $(subst ON,,$(MESG)),,$(shell mesr "$@\n$(GOAL_HERE)\nJUŻ\nKONIEC\nKOMPILACJI\n")))

# https://stackoverflow.com/questions/44488712/makefile-run-the-same-command-with-different-arguments-on-different-targets/44488843
# https://stackoverflow.com/questions/50322607/multiple-target-specific-variable-values

show_paths print_paths:
	@echo SOURCE_DIR       = $(SOURCE_DIR)
	@echo BUILD_PARENT_DIR = $(BUILD_PARENT_DIR)
	@echo BUILD_WORK_DIR   = $(BUILD_WORK_DIR)
	@echo INSTALL_DIR      = $(INSTALL_DIR)
	@echo THE_VERSION      = $(THE_VERSION)

# set variables for each prepare_* target
prepare_debug:        USE_DEBUG = 1
prepare_debug:        USE_FAST  = OFF
prepare_debug:        MAX_LOG   = 6
prepare_nodebug:      USE_DEBUG = 0
prepare_nodebug:      USE_FAST  = OFF
prepare_nodebug:      MAX_LOG   = 6
prepare_nodebug_fast: USE_DEBUG = 0
prepare_nodebug_fast: USE_FAST  = ON
prepare_nodebug_fast: MAX_LOG   = 3

# set all other cmake arguments
# remove Remove -DCHUNKSIZE=1 so that orignial file names are used.
#OTHER_CMAKE_ARGUMENTS    := -DENABLE_LOGGER=OFF

# Test all combinations: for ((i=1;i<19;i=i+1));do ; reset && echo FETS=$i && make uninstall && make FETS=$i && ./examples/yade --test && ./examples/yade --check || break ; done
#                        for ((i=1;i<19;i=i+1));do ; reset && echo FETS=$i && make uninstall && make FETS=$i && ./examples/yade --test && ./examples/yade --check || break ; done
# set -j number of jobs
JOBSNUM    := 32
VERBOSE    := 0
BITS       :=
DECI       := 15
MPFR       := 0
FETS       := 18
CRASH      := 0
SSE        := 0
THE_FAST   :=

ifeq ($(FETS),-15) # No DENABLE_POTENTIAL_BLOCKS
OTHER_CMAKE_ARGUMENTS     := -DENABLE_POTENTIAL_BLOCKS=OFF
endif
ifeq ($(FETS),-14) # ASAN
OTHER_CMAKE_ARGUMENTS     := -DENABLE_ASAN=1 -DENABLE_OPENMP=0
endif
ifeq ($(FETS),-13) # ASAN
OTHER_CMAKE_ARGUMENTS     := -DENABLE_ASAN=1
endif
ifeq ($(FETS),-12) # opposite
OTHER_CMAKE_ARGUMENTS     := -DENABLE_LIQMIGRATION=ON -DENABLE_PROFILING=ON -DENABLE_SPH=ON -DENABLE_DEFORM=ON
endif
ifeq ($(FETS),-11) # No PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_LINSOLV=OFF
endif
ifeq ($(FETS),-10) # No PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-10) # No PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-9) # No PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_GUI=OFF -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-8) # No PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_POTENTIAL_BLOCKS=OFF  -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GUI=OFF -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-7) # No VTK No GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_GUI=OFF -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-6) # No VTK:
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF
endif
ifeq ($(FETS),-5) # No VTK NO PP PB GUI No CGAL
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_POTENTIAL_BLOCKS=OFF  -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GUI=OFF -DENABLE_CGAL=OFF
endif
ifeq ($(FETS),-4) # No VTK NO PP PB GUI
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_POTENTIAL_BLOCKS=OFF  -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GUI=OFF
endif
ifeq ($(FETS),-3) # No VTK No PP PB
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_POTENTIAL_BLOCKS=OFF  -DENABLE_POTENTIAL_PARTICLES=OFF
endif
ifeq ($(FETS),-2) # Very Minimal + LOG ( 0):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_LOGGER=ON -DENABLE_OPENMP=OFF -DENABLE_USEFUL_ERRORS=0 -DENABLE_ASAN=1
endif
ifeq ($(FETS),-1) # Very Minimal + LOG ( 0):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_MPI=OFF -DENABLE_GUI=OFF -DENABLE_GTS=OFF -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=OFF -DENABLE_FEMLIKE=OFF -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=OFF -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GL2PS=OFF -DENABLE_LOGGER=ON -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),0) # Very Minimal + GUI + LOG ( 0):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_MPI=OFF -DENABLE_GUI=ON -DENABLE_GTS=OFF -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=OFF -DENABLE_FEMLIKE=OFF -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=OFF -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GL2PS=OFF -DENABLE_LOGGER=ON -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),1) # Very Minimal   ( 1):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=OFF -DENABLE_MPI=OFF -DENABLE_GUI=OFF -DENABLE_GTS=OFF -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=OFF -DENABLE_FEMLIKE=OFF -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=OFF -DENABLE_POTENTIAL_PARTICLES=OFF -DENABLE_GL2PS=OFF -DENABLE_LOGGER=OFF -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),2) # Middle         ( 2):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=ON  -DENABLE_MPI=OFF -DENABLE_GUI=OFF -DENABLE_GTS=ON  -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=OFF -DENABLE_FEMLIKE=ON  -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=ON  -DENABLE_POTENTIAL_PARTICLES=ON  -DENABLE_GL2PS=OFF -DENABLE_LOGGER=ON  -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),3) # Middle         ( 3):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=ON  -DENABLE_MPI=ON  -DENABLE_GUI=OFF -DENABLE_GTS=ON  -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=OFF -DENABLE_FEMLIKE=ON  -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=ON  -DENABLE_POTENTIAL_PARTICLES=ON  -DENABLE_GL2PS=ON  -DENABLE_LOGGER=ON  -DENABLE_OPENMP=ON
endif
ifeq ($(FETS),4) # CGAL           ( 4):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=ON  -DENABLE_MPI=ON  -DENABLE_GUI=OFF -DENABLE_GTS=ON  -DENABLE_TWOPHASEFLOW=OFF -DENABLE_CGAL=ON  -DENABLE_FEMLIKE=ON  -DENABLE_LBMFLOW=OFF -DENABLE_POTENTIAL_BLOCKS=ON  -DENABLE_POTENTIAL_PARTICLES=ON  -DENABLE_GL2PS=ON  -DENABLE_LOGGER=ON  -DENABLE_OPENMP=ON
endif
ifeq ($(FETS),5) # CGAL           ( 5):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_VTK=ON  -DENABLE_MPI=ON  -DENABLE_GUI=ON  -DENABLE_GTS=ON  -DENABLE_TWOPHASEFLOW=ON  -DENABLE_CGAL=ON  -DENABLE_FEMLIKE=ON  -DENABLE_LBMFLOW=ON  -DENABLE_POTENTIAL_BLOCKS=ON  -DENABLE_POTENTIAL_PARTICLES=ON  -DENABLE_GL2PS=ON  -DENABLE_LOGGER=ON  -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),6) # No OpenMP      ( 6):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_OPENMP=OFF
endif
ifeq ($(FETS),7) # No GUI         ( 7):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_GUI=OFF
endif
ifeq ($(FETS),8) # Mask Arbitrary ( 8):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_MASK_ARBITRARY=ON
endif
ifeq ($(FETS),9) # No Logger      ( 9):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_LOGGER=OFF
endif
# CHOLMOD_GPU SPH DEFORM LIQMIGRATION MASK_ARBITRARY PROFILING
ifeq ($(FETS),10) # SPH On        (10):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_SPH=ON
endif
ifeq ($(FETS),11) # DEFORM On     (11):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_DEFORM=ON
endif
ifeq ($(FETS),12) # LIQMIGRATION  (12):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_LIQMIGRATION=ON
endif
ifeq ($(FETS),13) # PROFILING On  (13):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_PROFILING=ON
endif
ifeq ($(FETS),14) # Defaults      (14):
OTHER_CMAKE_ARGUMENTS     := -DUSE_QT5=OFF
endif
ifeq ($(FETS),15) # Defaults      (15):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_MPI=OFF
endif
ifeq ($(FETS),16) # Defaults      (16):
OTHER_CMAKE_ARGUMENTS     := -DENABLE_MPI=OFF -DENABLE_MASK_ARBITRARY=ON
endif
ifeq ($(FETS),17) # Defaults      (17):
OTHER_CMAKE_ARGUMENTS     := -DDEBUG=ON
endif
ifeq ($(FETS),18) # Defaults      (18):
OTHER_CMAKE_ARGUMENTS     := #-DCMAKE_VERBOSE_MAKEFILE=ON
endif

## clang-8 : (zainstalowałem z debian buster)
#	export CC="clang-8" && \
#	export CXX="clang-8" && \
#	………… -DENABLE_USEFUL_ERRORS=0 …………
#


# those prepare_* targets all do the same thing, except the USE_DEBUG and USE_FAST are set differently
prepare_debug prepare_nodebug prepare_nodebug_fast:
	mkdir -p $(BUILD_PARENT_DIR) && \
	cd       $(BUILD_PARENT_DIR) && \
	mkdir -p $(BUILD_WORK_DIR) && \
	cd       $(BUILD_WORK_DIR) && \
	export CC="icc" && \
	export CXX="icpc" && \
	cmake -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR) $(SOURCE_DIR) \
	-DDEBUG=$(USE_DEBUG) -D$(THE_FAST) -DMAX_LOG_LEVEL=$(MAX_LOG) $(OTHER_CMAKE_ARGUMENTS) \
	-DCMAKE_VERBOSE_MAKEFILE=${VERBOSE} -DREAL_PRECISION_BITS=${BITS} -DREAL_DECIMAL_PLACES=${DECI} \
	-DENABLE_MPFR=${MPFR} -DALLOW_CRASHES=${CRASH} -DVECTORIZE=${SSE} -DENABLE_USEFUL_ERRORS=0

# -DEXPERIMENTALLY_FAST
# -DDANGEROUSLY_FAST
# declare an alias for "prepare" to invoke "prepare_nodebug"
prepare: prepare_nodebug

# those compilation targets all do the same thing, except the $@ evaluates to "all_debug", "all_nodebug", "all_nodebug_fast"
# and $(subst all_,,$@) removes the string "all_" from it, ( https://ftp.gnu.org/old-gnu/Manuals/make-3.79.1/html_chapter/make_8.html )
# and this is how the correct prepare_* target is called. But only if dir $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) does not exist.
# In case of error use echo to print a red+yellow message if there was an error. ( https://misc.flogisoft.com/bash/tip_colors_and_formatting )
all_debug all_nodebug all_nodebug_fast:
	test -d $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) || ${MAKE} prepare_$(subst all_,,$@)
	cd      $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) && \
	${MAKE} install -j $(JOBSNUM) \
	|| { echo "\e[91mERROR. \e[93mIf it's immediate try calling 'make prepare', or 'make clean'\e[0m, otherwise look for compilation errors." ; exit 1; }
	${MAKE} links
#	$(MAYBE_MESSAGE)

# declare an alias for "all" to invoke "all_nodebug"
all: all_nodebug
#	$(MAYBE_MESSAGE)

install: all_nodebug
#	$(MAYBE_MESSAGE)

docs_debug docs_nodebug docs_nodebug_fast:
	${MAKE} all_$(subst docs_,,$@)
	cd $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) && ${MAKE} doc -j $(JOBSNUM)
#	$(MAYBE_MESSAGE)

# declare an alias for "docs" to invoke "docs_nodebug"
docs: docs_nodebug
#	$(MAYBE_MESSAGE)

docsfast: all
	cd $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) && ${MAKE} doc/fast -j $(JOBSNUM)

run:
	$(INSTALL_DIR)/bin/yade-$(THE_VERSION)

clean:
	rm -rf $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR)

uninstall: clean
	rm -rf $(INSTALL_DIR)

# szybko odpalę terminale w tych katalogach. Przy okazji to sprawdza czy to są właściwe katalogi.
xbinstall:
	cd $(INSTALL_DIR) && xbterm

xbbuild:
	cd $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR) && xbterm

xbuild: xbbuild

xbmc:
	xbtermcomm mc ./ $(INSTALL_DIR)

# I use rox file manager for "clicking on files" :) You might want to replace this with `thunar` or something else.
# This command is great for debugging documentation.
roxdoc:
	rox $(INSTALL_DIR)/share/doc/yade-$(THE_VERSION)/html

roxbuild:
	rox $(BUILD_PARENT_DIR)/$(BUILD_WORK_DIR)

roxlib:
	rox $(INSTALL_DIR)/lib/x86_64-linux-gnu/yade-$(THE_VERSION)/

links:
	rm -f ./examples/yade-batch
	rm -f ./examples/yade
	ln -s $(INSTALL_DIR)/bin/yade-$(THE_VERSION)-batch ./examples/yade-batch
	ln -s $(INSTALL_DIR)/bin/yade-$(THE_VERSION)       ./examples/yade
	@echo THE_VERSION: $(THE_VERSION)

lsbin:
	ls -la $(INSTALL_DIR)/bin --color
	ls -la $(INSTALL_DIR)/lib/x86_64-linux-gnu --color
	ls -la ./examples/yade* --color


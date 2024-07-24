# ------------------------------------------------------------
# soft link  SCLCP src files to LAMMPS_SRC directory
# ------------------------------------------------------------

# can set LAMMPS_SRC from command line using "make LAMMPS_SRC=MYPATH"

LAMMPS_SRC ?= '../../../lammps_stable_16Mar2018/src/'

files = angle_orient.cpp \
		    angle_orient.h \
        angle_wlc_twist.cpp \
        angle_wlc_twist.h \
				fix_ext_alignment.cpp \
        fix_ext_alignment.h \

lammps_files := $(addprefix $(LAMMPS_SRC)/,$(files))
local_files := $(addprefix $(PWD)/,$(files))

# these are files that already exist in LAMMPS but must be overwritten in order for the model to compile
# typically they contain only minor changes from the original versions in order to suppress errors, etc
# Ideally there would be no files in this list, but in order to get things to work, we sometimes needed to hack LAMMPS a bit
override_file = atom.cpp

lammps_override_file := $(addprefix $(LAMMPS_SRC)/,$(override_file))
lammps_override_file_bak := $(addsuffix "_original",$(lammps_override_file))
local_override_file := $(addprefix $(PWD)/,$(override_file))


# link all .cpp .h files
#link:
#	ln -sfr *.cpp *.h -t ${LAMMPS_SRC}
#
#unlink:
#	find ${LAMMPS_SRC} -type l -delete


help:
	@echo "Please use one of the following options:"
	@echo "  make link       to soft link the model C++ files into the LAMMPS src directory "
	@echo "  make unlink     to remove the soft links between the model and LAMMPS "
	@echo ""
	@echo "Note: Path to LAMMPS source is  LAMMPS_SRC=${LAMMPS_SRC}"
	@echo "This can be changed using 'make link LAMMPS_SRC=<new path>.'"

# link only some files
link:
	@#ln -sfr $(files) -t ${LAMMPS_SRC} # Older ln doesnt permit the -r flag
	ln -sf $(local_files) -t ${LAMMPS_SRC}

	@#custom logic for atom.cpp (should remove eventually)
	@if [ `diff ${local_override_file} ${lammps_override_file} | grep "^>" | wc -l` -gt 2 ]; then \
		echo "${local_override_file} differs by more than two lines."; \
		exit 1; \
	else\
		echo "Note: Overriding ${lammps_override_file}. Saving old version to ${lammps_override_file_bak}";\
		cp ${lammps_override_file} ${lammps_override_file_bak};\
		ln -sf $(local_override_file) -t ${LAMMPS_SRC};\
	fi


unlink:
	rm ${lammps_files}

	@#restore original files that were overwritten
	mv ${lammps_override_file_bak} ${lammps_override_file}

#
# Makefile to compile the Fortran 90 source files using f2py
#

debug := 

RM := -rm
bold=`tput bold`
normal=`tput sgr0`

f90files := integratef90.f90

ifdef debug
    flags = "${my_gfortran_debug}"
else
    flags = ""
endif

all: $(f90files:.f90=.so)

%.so: %.f90
	@echo
	@echo "${bold}Building $@ ...${normal}"
	@f2py --fcompiler=gnu95 --f90flags=$(flags) -c $< -m $(basename $@) >> /dev/null
	@echo "${bold} ---Success---${normal}"

clean:
	$(RM) -f $(f90files:.f90=.so)
#$(RM) -f ./*.pyc

realclean: clean


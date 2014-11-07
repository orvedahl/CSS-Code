#
# Makefile to compile the Fortran 90 source files using f2py (and gfortran)
#

debug := 

RM := -rm
bold=`tput bold`
normal=`tput sgr0`

f90files := integratef90.f90 shell_avgf90.f90

ifdef debug
    flags = "${my_gfortran_debug}"
else
    flags = ""
endif

all: $(f90files:.f90=.so)

# explicitly put in the individual dependencies (kludgy, but it works)
# ORDER MATTERS (because these are dependencies)
integratef90.so: integratef90.f90
shell_avgf90.so: integratef90.f90 shell_avgf90.f90

# the "$^" macro holds the dependencies (and their full path)
%.so:
	@echo
	@echo "${bold}Building $@ ...${normal}"
	@f2py --fcompiler=gnu95 --f90flags=$(flags) -c $^ -m $(basename $@) >> /dev/null
	@echo "${bold} ---Success---${normal}"

clean:
	$(RM) -f $(f90files:.f90=.so)

realclean: clean

# debugging tool
print-%: ; @echo $* is $($*)

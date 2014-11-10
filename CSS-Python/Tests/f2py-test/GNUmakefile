#
# Makefile to compile the Fortran 90 source files & compile using F2py
#

python := t

debug := 

##########################################################################
# Fortran source files, and which files to build using f2py
f90files := integratef90.f90 functionf90.f90 fortran_driver.f90
pyfiles := integratef90.f90 functionf90.f90
##########################################################################

RM := -rm
bold=`tput bold`
normal=`tput sgr0`

# compiler flags
ifdef debug
    flags = ${my_gfortran_debug}
    f2pyflags = "${my_gfortran_debug}"  # f2py requires flags as a string
else
    flags = -O2 -fno-range-check
    f2pyflags = ""
endif

ifdef python
all: xtest.exe $(pyfiles:.f90=.so)
else
all: xtest.exe
endif

##########################################################################
# explicitly put in the individual dependencies (kludgy, but it works)
# ORDER MATTERS!!
xtest.exe: integratef90.o functionf90.o fortran_driver.o
integratef90.o: integratef90.f90
functionf90.o: integratef90.o functionf90.f90
fortran_driver.o: functionf90.o fortran_driver.f90

# include dependencies for f2py
ifdef python
integratef90.so: integratef90.f90
functionf90.so: integratef90.f90 functionf90.f90
endif
##########################################################################

xtest.exe: $(f90files:.f90=.o)
	@echo
	@echo "${bold}Building $@ ...${normal}"
	gfortran $(flags) $^ -o $@
	@echo "${bold} ---Success---${normal}"

%.o: %.f90
	@echo
	@echo "${bold}Building $@ ...${normal}"
	gfortran $(flags) -c $< -o $@
	@echo "${bold} ---Success---${normal}"

%.so:
	@echo
	@echo "${bold}Building $@ ...${normal}"
	@f2py --fcompiler=gnu95 --f90flags=$(f2pyflags) -c $^ -m $(basename $@) >> /dev/null
	@echo "${bold} ---Success---${normal}"

clean:
	$(RM) -f $(f90files:.f90=.o)
	$(RM) -f $(f90files:.f90=.mod)

realclean: clean
	$(RM) -f xtest.exe
	$(RM) -f $(f90files:.f90=.so)

# debugging tool
print-%: ; @echo $* is $($*)

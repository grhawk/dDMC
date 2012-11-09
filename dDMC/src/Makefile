PROG = dDMC
FC = ifort
#FC = gfortran
#FC = gfortran -ffixed-line-length-0 -std=legacy 
LIBS = 
FFLAGS = -O -g
LDFLAGS = -g #-s

SRCS =  Makefile file_tools.f90 precision.f90 main.f90 read_atomdata.f90 \
	read_tag_dftbp.f90 read_xyz.f90 string_tools.f90 utils.f90 parameters.f90\
	ReadInput.f90

OBJS =  main.o file_tools.o precision.o read_atomdata.o \
	read_tag_dftbp.o read_xyz.o string_tools.o utils.o parameters.o \
	ReadInput.o

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

# use this PATTERN rule instead of old-fashioned SUFFIX rules like 
#.SUFFIXES: $(SUFFIXES) .f90 
#
#.f90.o:
#	$(FC) $(FFLAGS) -c $<
%.o : %.f90
	$(FC) $(FFLAGS) -c $<


file_tools.o: precision.o
main.o: precision.o read_atomdata.o read_tag_dftbp.o read_xyz.o utils.o parameters.o ReadInput.o file_tools.o
read_atomdata.o: file_tools.o precision.o
read_tag_dftbp.o: file_tools.o precision.o string_tools.o
read_xyz.o: file_tools.o precision.o
string_tools.o: precision.o
parameters.o: precision.o
utils.o: precision.o
ReadInput.o: precision.o string_tools.o

.PHONY: default clean veryclean debug fcheck fdepend showtargets

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod $(PROG)

# veryclean: remove all regenerable files, including executables
veryclean: clean
	$(RM) $(PROG) $(OBJS)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g " "OSTYPE=$(OSTYPE:-gnu=)"

# fcheck: pass sources to the Fortran program checker
fcheck:
	fcheck -nonovice main.f90

# fdepend: pass sources to the Fortran dependency analyzer
fdepend:
	fdepend main.f90

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean veryclean debug fcheck fdepend 

distclean: veryclean

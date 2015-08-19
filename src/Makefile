PROG = dDMC2
#FC = ifort
FC = gfortran
LIBS = 
FFLAGS = -O
LDFLAGS = 

SRCS =  Makefile file_tools.f90 precision.f90 main.f90 read_atomdata.f90 \
	read_tag_dftbp.f90 read_xyz.f90 string_tools.f90 utils.f90 parameters.f90\
	read_input.f90 dampingfunctions.f90 read_gf.f90

OBJS =  main.o file_tools.o precision.o read_atomdata.o \
	read_tag_dftbp.o read_xyz.o string_tools.o utils.o parameters.o \
	read_input.o dampingfunctions.o read_gf.o

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

file_tools.o: precision.o
main.o: precision.o read_atomdata.o read_tag_dftbp.o read_xyz.o utils.o \
	parameters.o read_input.o file_tools.o dampingfunctions.o \
	read_gf.o Makefile
read_atomdata.o: file_tools.o precision.o
read_tag_dftbp.o: file_tools.o precision.o string_tools.o
read_xyz.o: file_tools.o precision.o
string_tools.o: precision.o
parameters.o: precision.o
utils.o: precision.o
read_input.o: precision.o string_tools.o
dampingfunctions.o: precision.o file_tools.o read_input.o

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

fullclean: veryclean

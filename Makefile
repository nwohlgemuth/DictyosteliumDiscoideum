default: fluid

help:
	@echo 'Try one of:'
	@echo '  make cellmove      - make the actual program'
	@echo '  make clean      - clean up all files from compiling & running'
	@echo '  make linecount  - write the number of lines in the program'
	@echo '  make sloccount  - statistics on source lines of code count'

.PHONY : clean linecount testing

# BASIC_WARNINGS -- keeps the compile line simple
# MAIN_WARNINGS  -- all warnings I can expect my code to not ignore
# MEGA_WARNINGS  -- extra warnings that can help, but might cause spuriousness
# WARNINGS       -- warnings level I actually use
# OPTIMIZATIONS  -- optimization flags to compile with
# PROFILING      -- profiling flag to compile with (so gprof can be used)
# TESTING_FLAGS  -- testing flags to use (don't set manually, see testing target)
BASIC_WARNINGS = -Wall -Werror
MAIN_WARNINGS = -Wall -W -Wabi -Wcast-align -Wcast-qual -Wconversion \
                -Wdisabled-optimization -Wendif-labels -Winline \
                -Wpointer-arith -Wredundant-decls -Wshadow \
                -Wsign-promo -Wundef -Wwrite-strings -Werror
MEGA_WARNINGS = $(MAIN_WARNINGS) \
                -Waggregate-return -Weffc++ -Wpadded -Wunreachable-code \
	        -Wold-style-cast -Woverloaded-virtual -Wmissing-noreturn 
WARNINGS      = $(BASIC_WARNINGS)
OPTIMIZATIONS = -O3 -funroll-loops # -fno-inline -Wno-inline ## super-debugging
PROFILING = -p
LOG_FLAGS = -DI_WANT_SPEED  # defining I_WANT_SPEED turns logging off
nTESTING_FLAGS = # This should be left blank; it's set manually for testing
SLEPC_INCLUDES=-I/usr/bmake -I/usr/include/mpiuni 
SLEPC_LIBS=-lslepc -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat \
           -lpetscvec -lpetsc -lmpiuni -llapack -lblas -lm
ARPACK_LIBS=-larpack
LDFLAGS = #-framework vecLib
	# Leaving this blank for now; may add -rdynamic for testing; the
          # -rdynamic flag is needed for util.cc:show_stack_trace to provide
          # useful output
LDLIBS  = #-L/sw/lib -lfftw3
	  # $(ARPACK_LIBS) $(SLEPC_LIBS)
          # -llapack

#CXXFLAGS= -pipe $(WARNINGS) $(TESTING_FLAGS) $(SLEPC_INCLUDES) \
#	  $(OPTIMIZATIONS) $(LOG_FLAGS) # $(PROFILING)
CXXFLAGS= $(OPTIMIZATIONS)

OBJS = cellclass.o usecell.o mersenne.o 
OBJNS =  mersenne.o
CFILES = cellclass_simple_1.cc stiffusecell.cc stoc1.cpp stoc2.cpp userintf.cpp cellclass.h decvar.h


cellmoveclang: $(CFILES)
	g++ cellclass_simple_1.cc stiffusecell.cc stoc1.cpp stoc2.cpp userintf.cpp  $(CXXFLAGS) $(LDFLAGS) $(OBJNS) -stdlib=libc++ -std=c++11 -lc++ -lsundials_cvode -lsundials_nvecserial -o $@

clean:
	$(RM) fluid testing *.o *~

linecount:
	wc -l *.I *.h *.cc

sloccount:
	sloccount .  # --details before . is allowed

filecount:
	@# Notes:
	@# (1) In Makefiles, I'm forced to use '$$' when I mean '$' in order to
	@#     tell the make program that I really do want a $ passed to bash;
	@#     otherwise make stupidly tries to find some weird variable and
	@#     errors out on me.
	@# (2) By having make invoke itself, there are two extra lines of
	@#     output; combined with the wc command and its that adds up to 4
	@#     lines that we have to subtract off.
	@echo $$(( $$(make linecount | wc -l) - 4 ))
	@# Easier way (but who wants to do things the easy way?):
	@#   ls *.I *.h *.cc | wc -l

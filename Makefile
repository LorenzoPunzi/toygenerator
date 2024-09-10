# Makefile for compiling .f90 files using gfortran

# Compiler and compiler flags
CXX = gfortran

# Source files (all .f90 files in the current directory)
SRCS = $(wildcard *.f90)

# Executable names (remove .f90 extension)
EXECS = $(patsubst %.f90,%,$(SRCS))

all: $(EXECS)

%: %.f90
	$(CXX) -o $@ $<

clean:
	rm -f $(EXECS)

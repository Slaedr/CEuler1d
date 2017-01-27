# Makefile for Euler1d project
# Aditya Kashi
#
# NOTE: Make sure environment variable CC has been set with the compilers to use, or set them below.
# NOTE: -Msafeptr=all means all pointers are assumed restriced.

NAME = euler1d-acc
PREFIX = build

INCLUDES = 

#Profile using gprof
PROFILE= #-pg

CFLAGS =

CLIBS = -lm

# if DEBUG is 1 or not defined, code is compiled in debug mode. Otherwise, optimizations are enabled
ifndef DEBUG

  $(info "Compiling with optimizations, without debug data")
  ifeq ($(CC),pgcc)
    $(info "Setting flags for pgc++")
    CFLAGS = -O3 -Msafeptr=all -fast #-Minfo=vect -Minfo=inline
    LFLAGS = -O3
  else
    CFLAGS =  -O3 -Winline -ftree-vectorizer-verbose=2
    LFLAGS = -O3
  endif

else

  #PROFILE = -pg
  $(info "Compiling debug version")
  ifeq ($(CC),pgcc)
    CFLAGS = -g #-Minfo=vect,inline
    LFLAGS = 
  else
    CFLAGS = -ggdb -Winline -ftree-vectorizer-verbose=2
    LFLAGS = -ggdb
  endif

endif

ifdef BUILD_WITH_ACC
  $(info 'Compiling with OpenACC')
  ifeq ($(CXX),pgcc)
    ifdef BUILD_FOR_MULTICORE
      $(info 'Compiling for multicore CPU')
      CFLAGS := $(CFLAGS) -ta=multicore -Minfo=accel
      LFLAGS := $(LFLAGS) -ta=multicore
    else
      $(info 'Compiling for Nvidia GPU')
      CFLAGS := $(CFLAGS) -ta=tesla -Minfo=accel
      LFLAGS := $(LFLAGS) -ta=tesla
    endif
  endif
endif


clibsrcs =$(wildcard *.c)
clibobjst =$(clibsrcs:.c=.o)
clibobjs = $(foreach obj,$(clibobjst),$(PREFIX)/$(obj))

$(NAME): $(clibobjs)
	$(CC) $(LFLAGS) -o $(PREFIX)/$(NAME) $(clibobjs) $(CLIBS) $(PROFILE)

$(PREFIX)/%.o: %.c
	$(CC)  $(CFLAGS) -c -o $@ $<  $(INCLUDES) $(CLIBS) $(PROFILE)

.PHONY : clean
clean:
	rm -f $(PREFIX)/$(NAME) $(libobjs) $(clibobjs)
	rm -f ./output/*


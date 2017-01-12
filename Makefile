# Makefile for Euler1d project
# Aditya Kashi
#
# Make sure environment variables CC and CXX have been set with the compilers to use, or set them below.

NAME = euler1d-acc
PREFIX = build

INCLUDES = 

#Profile using gprof
PROFILE= #-pg

CFLAGS =

# if DEBUG is 1 or not defined, code is compiled in debug mode. Otherwise, optimizations are enabled
ifndef DEBUG

  $(info "Compiling with optimizations, without debug data")
  ifeq ($(CXX),pgc++)
    $(info "Setting flags for pgc++")
    CPPFLAGS = -std=c++11 -O3 -Msafeptr=all -fast #-Minfo=vect -Minfo=inline
    LFLAGS = -O3
  else
    CPPFLAGS =  -std=c++14 -O3 -Winline -ftree-vectorizer-verbose=2
    LFLAGS = -O3 #-lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lpthread
  endif

else

  PROFILE = -pg
  $(info "Compiling debug version")
  ifeq ($(CXX),pgc++)
    CPPFLAGS = -std=c++11 -g #-Minfo=vect,inline
    LFLAGS = 
  else
    CPPFLAGS =  -std=c++14 -ggdb -Winline -ftree-vectorizer-verbose=2
    LFLAGS = -ggdb #-lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lpthread
  endif

endif

ifdef BUILD_WITH_ACC
  $(info 'Compiling with OpenACC')
  ifeq ($(CXX),pgc++)
    ifdef BUILD_FOR_MULTICORE
      $(info 'Compiling for multicore CPU')
      CPPFLAGS := $(CPPFLAGS) -ta=multicore -Minfo=accel
      LFLAGS := $(LFLAGS) -ta=multicore
    else
      $(info 'Compiling for Nvidia GPU')
      CPPFLAGS := $(CPPFLAGS) -ta=tesla -Minfo=accel
      LFLAGS := $(LFLAGS) -ta=tesla
    endif
  endif
endif
 
libsrcs =$(wildcard *.cpp)     
libobjst = $(libsrcs:.cpp=.o)
libobjs = $(foreach obj,$(libobjst),$(PREFIX)/$(obj))
clibsrcs =$(wildcard *.c)
clibobjst =$(clibsrcs:.c=.o)
clibobjs = $(foreach obj,$(clibobjst),$(PREFIX)/$(obj))

$(NAME): $(libobjs)
	$(CXX) $(LFLAGS) -o $(PREFIX)/$(NAME) $(libobjs) $(LIBS) $(PROFILE)

$(PREFIX)/%.o: %.cpp                    
	$(CXX) $(CPPFLAGS) $(EIGENDEF) -c -o $@ $<  $(INCLUDES) $(PROFILE)

$(PREFIX)/%.o: %.c
	$(CC)  $(CFLAGS) -c -o $@ $<  $(INCLUDES) $(PROFILE)

.PHONY : clean
clean:
	rm -f $(PREFIX)/$(NAME) $(libobjs) $(clibobjs)
	rm -f ./output/*


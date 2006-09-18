MAINTAINERCLEANFILES = Makefile.in

if WANT_GNUSEARCH_H
  HSEARCH_R_INC_DIR = -I$(top_srcdir)/lib
else
  HSEARCH_R_INC_DIR = 
endif

INCLUDES = -I$(top_builddir) -I$(top_srcdir) $(HSEARCH_R_INC_DIR)

lib_LIBRARIES       = libSNet.a 
libSNet_a_SOURCES = barrier.c barrier.h element.cc element.h errors.h nlayer.cc nlayer.h nnet.cc nnet.h progobj.cc progobj.h socketobj.cc socketobj.h timers.cc timers.h

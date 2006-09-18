MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = no-dependencies

if WANT_HSEARCH_R_O
  HSEARCH_R_OBJ = libhsearch_r.a
else
  HSEARCH_R_OBJ =
endif

noinst_LIBRARIES = $(HSEARCH_R_OBJ) libgetopt.a

libhsearch_r_a_SOURCES = hsearch_r.cpp gnusearch.h
libhsearch_r_a_CXXFLAGS = -I. -I$(top_srcdir) -I$(top_builddir)

libgetopt_a_SOURCES = getopt.c getopt.h
libgetopt_a_CFLAGS = -I$(top_srcdir) -I$(top_builddir)



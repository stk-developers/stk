
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = STKLib SNetLib .

if WANT_GNUSEARCH_H
  HSEARCH_R_INC_DIR = -I$(top_srcdir)/lib
else
  HSEARCH_R_INC_DIR =
endif

if WANT_LIBGETOPT
  LIBGETOPT = -lgetopt
else
  LIBGETOPT =
endif

# These are the common libraries needed by all programs
COMMONLDFLAGS = -pthread
COMMONLIBS    = ./STKLib/libSTKLib.a $(LIBGETOPT)

AM_LDFLAGS  = $(COMMONLDFLAGS) 
LDADD       = $(COMMONLIBS) 


################################################################################
INCLUDES    = $(ATLASINCLUDEDIR) -I. -I$(top_builddir) -I$(top_srcdir) \
              $(HSEARCH_R_INC_DIR)


# These are the programs that we want to build in this directory
bin_PROGRAMS = SExpand SERest SVite SNet SFeaCat SLRatio SResults SRover


# Specifications
SERest_SOURCES    = SERest.cc 
SVite_SOURCES     = SVite.cc
SExpand_SOURCES   = SExpand.cc
SFeaCat_SOURCES   = SFeaCat.cc
SLRatio_SOURCES   = SLRatio.cc
SResults_SOURCES  = SResults.cc
SRover_SOURCES    = SRover.cc

# SNet is kind of special. It needs some files in the SNet source. This will probably change and separate SNetLib will be built
SNet_SOURCES      = SNet.cc 
SNet_LDFLAGS      = $(COMMONLDFLAGS) -pthread
SNet_LDADD        = $(COMMONLIBS) ./SNetLib/libSNet.a

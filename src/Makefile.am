
MAINTAINERCLEANFILES = Makefile.in

if USE_SNET
  SNETLIBF = SNetLib
	SNETF    = SNet
else
  SNETLIBF = 
	SNETF    = 
endif

SUBDIRS = STKLib $(SNETLIBF) .

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

if WANT_HSEARCH_R_O
  LIBHSEARCH_R = -lhsearch_r
else 
  LIBHSEARCH_R = 
endif

if ENABLE_SSE
  SSE_FLAGS = -msse -msse2 -mmmx -mfpmath=sse
else
  SSE_FLAGS =
endif

AM_CXXFLAGS = $(CXXFLAGS) $(SSE_FLAGS)
AM_CFLAGS   = $(CFLAGS) $(SSE_FLAGS)

# These are the common libraries needed by all programs
COMMONLDFLAGS = -pthread -L$(top_builddir)/lib
COMMONLIBS    = ./STKLib/libSTKLib.a $(LIBGETOPT) $(LIBHSEARCH_R) 

AM_LDFLAGS  = $(COMMONLDFLAGS) 
LDADD       = $(COMMONLIBS) 


################################################################################
INCLUDES    = $(ATLASINCLUDEDIR) -I. -I$(top_builddir) -I$(top_srcdir) \
              $(HSEARCH_R_INC_DIR)


# These are the programs that we want to build in this directory
# bin_PROGRAMS = SExpand SERest SVite $(SNETF) SFeaCat SLRatio SResults SRover 
if USE_SNET
bin_PROGRAMS = SExpand SERest SVite SFeaCat SFeaStack SLRatio SResults SRover STree SNet SNetStats SFeaStat
else
bin_PROGRAMS = SExpand SERest SVite SFeaCat SFeaStack SLRatio SResults SRover STree SFeaStat
endif


# Specifications
SERest_SOURCES    = SERest.cc 
SVite_SOURCES     = SVite.cc
SExpand_SOURCES   = SExpand.cc
SFeaCat_SOURCES   = SFeaCat.cc
SFeaStack_SOURCES   = SFeaStack.cc
SLRatio_SOURCES   = SLRatio.cc
SResults_SOURCES  = SResults.cc
SRover_SOURCES    = SRover.cc
STree_SOURCES     = STree.cc
STree_LDADD       = $(COMMONLIBS) $(BOOST_PROGRAM_OPTIONS_LIB)
SFeaStat_SOURCES  = SFeaStat.cc
SFeaStat_LDADD    = $(COMMONLIBS) $(BOOST_PROGRAM_OPTIONS_LIB)
# SWEChCoFast_SOURCES = SWEChCoFast.cc
# SWE_SOURCES       = SWE.cc

#SERest_DEPENDENCIES    = $(RECURSIVE_TARGETS) 
#SVite_DEPENDENCIES     = $(RECURSIVE_TARGETS)
#SExpand_DEPENDENCIES   = $(RECURSIVE_TARGETS)
#SFeaCat_DEPENDENCIES   = $(RECURSIVE_TARGETS)
#SLRatio_DEPENDENCIES   = $(RECURSIVE_TARGETS)
#SResults_DEPENDENCIES  = $(RECURSIVE_TARGETS)
#SRover_DEPENDENCIES    = $(RECURSIVE_TARGETS)

# SNet is kind of special. It needs some files in the SNet source. This will probably change and separate SNetLib will be built
if USE_SNET
SNet_SOURCES      = SNet.cc 
SNet_LDFLAGS      = $(COMMONLDFLAGS) -pthread
SNet_LDADD        = $(COMMONLIBS) ./SNetLib/libSNet.a

SNetStats_SOURCES      = SNetStats.cc 
SNetStats_LDFLAGS      = $(COMMONLDFLAGS) -pthread
SNetStats_LDADD        = $(COMMONLIBS) ./SNetLib/libSNet.a
endif

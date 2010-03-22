MAINTAINERCLEANFILES = Makefile.in

if WANT_GNUSEARCH_H
  HSEARCH_R_INC_DIR = -idirafter $(top_srcdir)/lib
else
  HSEARCH_R_INC_DIR = 
endif


if WANT_HSEARCH_R_O
  HSEARCH_R_OBJ = $(top_srcdir)/lib/libhsearch_r.a
else 
  HSEARCH_R_OBJ = 
endif

if ENABLE_SSE
  SSE_FLAGS = -msse -msse2 -mmmx -mfpmath=sse
else
  SSE_FLAGS =
endif


INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
            $(HSEARCH_R_INC_DIR) @ATLASINCLUDEDIR@  \
						@MATLABINCLUDEDIR@

lib_LIBRARIES       = libSTKLib.a 
libSTKLib_a_LIBADD  = $(HSEARCH_R_OBJ)
libSTKLib_a_CXXFLAGS = $(SSE_FLAGS) $(AM_CXXFLAGS)
libSTKLib_a_CFLAGS   = $(SSE_FLAGS) $(AM_CFLAGS)
EXTRA_libSTKLib_a_SOURCES  = $(IMAGESC_SRC)
libSTKLib_a_DEPENDENCIES = $(HSEARCH_R_OBJ)
libSTKLib_a_SOURCES = BasicVector.cc BasicVector.h BasicVector.tcc \
		BDTree.h BDTree.cc BDTree.tcc BDTree_IO.cc BDTree_IO.h \
		BQuestion.h \
		ContextSample.cc ContextSample.h Error.cc  \
		Error.h Features.cc Features.h Matrix.cc Matrix.h Matrix.tcc \
		MlfStream.h MlfStream.cc MlfStream.tcc Models.cc  \
		Models.h Models_IO.cc Net.h Net.cc Net.tcc Decoder.cc Decoder.tcc \
		Decoder.h DecoderNetwork.cc DecoderNetwork.h DecoderNetwork_IO.cc DecoderNetwork_IO.tcc ZNorm.cc  \
		ZNorm.h common.C common.h dict.C dict.h fileio.h filmatch.C filmatch.h strtod_lc.C strtod_lc.h  \
		fileio.C getopt.C getopt.h labels.C labels.h stkstream.h stkstream.tcc   \
		Tokenizer.h Tokenizer.cc \
		Lattice.cc Lattice.h SigP.cc SigP.h \
		FeaCatAdapter.cc FeaCatAdapter.h

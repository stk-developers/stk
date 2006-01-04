###############################################################################
# compiler options
CFLAGS = -I/usr/X11R6/include -Wall -g -DDOUBLEPRECISION=0  
CXXFLAGS = $(CFLAGS)

###############################################################################
# linker options
LDLIBS = -lm -latlas -lpthread
#LDFLAGS = -a -g -pg

###############################################################################
# -DEXPANDNET_ONLY
XFLAGS = -DMOTIF -L/usr/X11R6/lib -lXm -lXt -lX11 -lXpm

#ICCFLAGS = -O2 -tpp6 -axiMK -xiMK -parallel

###############################################################################
# compiler binary
CC  = g++
CXX = g++

###############################################################################
# tools summary (for make all)
TOOLS =     SVite \
            SRover \
            SERest \
            SResults \
            SExpand \
            SLRatio


LIBHMMS = STKLib/Models.o STKLib/Matrix.o

###############################################################################
# object needed by all tools
BASELIBS =  STKLib/common.o \
            STKLib/Error.o \
            STKLib/filmatch.o
            
###############################################################################
# implicit "make all" rule             
all  : $(TOOLS)


###############################################################################
# tools-specific rules
SFea2OutP   : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/net.o
SVite       : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o
SLRatio     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o
SERest      : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o 
SFeaCat     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS)
SRover      : $(BASELIBS) STKLib/labels.o
SResults    : $(BASELIBS) STKLib/labels.o
SExpand     : $(BASELIBS) STKLib/labels.o STKLib/dict.o STKLib/net.o

SERestX :  SERest.o $(BASELIBS) $(LIBHMMS) STKLib/viterbi.c STKLib/labels.o STKLib/net.o STKLib/dict.o STKLib/imagesc.o
	$(CC) $(CFLAGS) $(XFLAGS) -o SERestX SERest.o STKLib/common.o STKLib/fileio.o STKLib/Models.o STKLib/viterbi.c STKLib/labels.o STKLib/net.o STKLib/dict.o STKLib/imagesc.o -lm


###############################################################################
# library specific rules
#


###############################################################################
# cleans the binaries  
clean:
	-rm -f *.o STKLib/*.o SVite SFea2OutP SRover SERest SLRatio SResults SExpand SCombAcc SMMIUpdate


###############################################################################
# runs makedepend to automatically define header files dependencies
depend:
	-makedepend -Y *.c *.cc STKLib/*.c STKLib/*.cc


###############################################################################
# some testing rules
test : $(BASELIBS)
test2 : test2.cc STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/Buffer.o

# DO NOT DELETE

SERest.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SERest.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SERest.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/viterbi.h
SERest.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/net.h
SERest.o: STKLib/labels.h STKLib/stkstream.h STKLib/MmfIO.h
SExpand.o: STKLib/net.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SExpand.o: STKLib/Error.h STKLib/dict.h STKLib/labels.h STKLib/common.h
SFeaCat.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFeaCat.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SFeaCat.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/viterbi.h
SFeaCat.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/net.h
SFea2OutP.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFea2OutP.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SFea2OutP.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/viterbi.h
SFea2OutP.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/net.h
SLRatio.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SLRatio.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/Models.h
SLRatio.o: STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SLRatio.o: STKLib/stkstream.tcc STKLib/fileio.h STKLib/labels.h
SLRatio.o: STKLib/common.h
SResults.o: STKLib/labels.h STKLib/common.h
SRover.o: STKLib/labels.h STKLib/common.h
SVite.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SVite.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/Models.h
SVite.o: STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SVite.o: STKLib/stkstream.tcc STKLib/fileio.h STKLib/labels.h STKLib/common.h
SNet.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SNet.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SNet.o: STKLib/stkstream.tcc STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
SNet.o: STKLib/dict.h STKLib/net.h
test.o: STKLib/stkstream.h
test2.o: STKLib/Matrix.h STKLib/Error.h STKLib/Matrix.tcc STKLib/common.h
test2.o: STKLib/String.h STKLib/Buffer.h STKLib/stkstream.h
STKLib/common.o: STKLib/common.h STKLib/Error.h STKLib/filmatch.h
STKLib/dict.o: STKLib/dict.h STKLib/Models.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/stkstream.h STKLib/stkstream.tcc
STKLib/filmatch.o: STKLib/filmatch.h STKLib/common.h STKLib/Error.h
STKLib/getopt.o: STKLib/getopt.h
STKLib/imagesc.o: STKLib/imagesc.h
STKLib/labels.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
STKLib/labels.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/net.h
STKLib/net.o: STKLib/Models.h STKLib/labels.h STKLib/common.h STKLib/Error.h
STKLib/net.o: STKLib/dict.h STKLib/net.h
STKLib/viterbi.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
STKLib/viterbi.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/net.h
STKLib/Buffer.o: STKLib/Buffer.h
STKLib/Error.o: STKLib/Error.h
STKLib/labelreader.o: STKLib/labelreader.h STKLib/stkstream.h STKLib/common.h
STKLib/labelreader.o: STKLib/Error.h STKLib/stkstream.tcc
STKLib/Matrix.o: STKLib/Matrix.h STKLib/Error.h STKLib/Matrix.tcc
STKLib/Matrix.o: STKLib/common.h
STKLib/MmfIO.o: STKLib/MmfIO.h STKLib/Models.h STKLib/Error.h
STKLib/Models.o: STKLib/Models.h STKLib/stkstream.h STKLib/common.h
STKLib/Models.o: STKLib/Error.h STKLib/stkstream.tcc

CFLAGS = -I/usr/X11R6/include -Wall -g -DDOUBLEPRECISION=0
# -DEXPANDNET_ONLY
XFLAGS = -DMOTIF -L/usr/X11R6/lib -lXm -lXt -lX11 -lXpm
#LDFLAGS = -a -g -pg
#ICCFLAGS = -O2 -tpp6 -axiMK -xiMK -parallel
CC  = g++


###############################################################################
# tools summary (for make all)
TOOLS =     SVite \
            SRover \
            SERest \
            SResults \
            SExpand \
            SLRatio


LIBHMMS = STKLib/Models.o

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
SFea2OutP   : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/net.o -lm
SVite       : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SLRatio     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SERest      : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SFeaCat     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) -lm
SRover      : $(BASELIBS) STKLib/labels.o -lm
SResults    : $(BASELIBS) STKLib/labels.o -lm
SExpand     : $(BASELIBS) STKLib/labels.o STKLib/dict.o STKLib/net.o -lm

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
	-makedepend -Y *.c HTKLib/*.c STKLib/*.c STKLib/*.cc


###############################################################################
# some testing rules
test : test.cc STKLib/common.o STKLib/labelreader.o STKLib/filmatch.o STKLib/Buffer.o
test2 : test2.cc STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/Buffer.o


# DO NOT DELETE

SERest.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SERest.o: STKLib/Models.h STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
SERest.o: STKLib/dict.h STKLib/net.h STKLib/labels.h
SExpand.o: STKLib/net.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SExpand.o: STKLib/Error.h STKLib/dict.h STKLib/labels.h STKLib/common.h
SFeaCat.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFeaCat.o: STKLib/Models.h STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
SFeaCat.o: STKLib/dict.h STKLib/net.h
SFea2OutP.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFea2OutP.o: STKLib/Models.h STKLib/viterbi.h STKLib/Models.h STKLib/labels.h
SFea2OutP.o: STKLib/dict.h STKLib/net.h
SLRatio.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SLRatio.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/Models.h
SLRatio.o: STKLib/fileio.h STKLib/labels.h STKLib/common.h
SResults.o: STKLib/labels.h STKLib/common.h
SRover.o: STKLib/labels.h STKLib/common.h
SVite.o: STKLib/viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SVite.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/Models.h
SVite.o: STKLib/fileio.h STKLib/labels.h STKLib/common.h
STKLib/common.o: STKLib/common.h STKLib/Error.h STKLib/filmatch.h
STKLib/dict.o: STKLib/dict.h STKLib/Models.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h
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
STKLib/labelreader.o: STKLib/labelreader.h STKLib/stkstream.h
STKLib/labelreader.o: STKLib/stkstream.tcc STKLib/common.h STKLib/Error.h
STKLib/Models.o: STKLib/Models.h STKLib/stkstream.h STKLib/stkstream.tcc
STKLib/Models.o: STKLib/common.h STKLib/Error.h
STKLib/ModelSet.o: STKLib/stkstream.h STKLib/stkstream.tcc

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

XXX = `ls`

###############################################################################
# tools summary (for make all)
TOOLS =     SVite \
            SRover \
            SERest \
            SResults \
            SExpand \
            SFeaCat \
            SLRatio \
	    SNet


###############################################################################
# library-specific modules
LIBHMMS = STKLib/Models.o STKLib/Models_IO.o STKLib/Matrix.o
LIBNET  = STKLib/Net.o STKLib/Net_IO.o


###############################################################################
# object needed by all tools
BASELIBS =  STKLib/common.o \
            STKLib/Error.o \
            STKLib/filmatch.o
            
###############################################################################
# implicit "make all" rule             
all  : $(TOOLS)
	-echo "Done :-)"

###############################################################################
# tools-specific rules
SFea2OutP   : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.o $(LIBNET)
SVite       : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.o $(LIBNET) STKLib/labels.o STKLib/dict.o
SLRatio     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.o $(LIBNET) STKLib/labels.o STKLib/dict.o
SERest      : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.o $(LIBNET) STKLib/labels.o STKLib/dict.o 
SFeaCat     : $(BASELIBS) STKLib/fileio.o $(LIBHMMS)
SRover      : $(BASELIBS) STKLib/labels.o
SResults    : $(BASELIBS) STKLib/labels.o
SExpand     : $(BASELIBS) STKLib/labels.o STKLib/dict.o $(LIBNET)

SERestX :  SERest.o $(BASELIBS) $(LIBHMMS) STKLib/Viterbi.cc STKLib/labels.o STKLib/Net.o STKLib/dict.o STKLib/imagesc.o
	$(CC) $(CFLAGS) $(XFLAGS) -o SERestX SERest.o STKLib/common.o STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.cc STKLib/labels.o $(LIBNET) STKLib/dict.o STKLib/imagesc.o -lm
	
SNet : $(BASELIBS) STKLib/fileio.o $(LIBHMMS) STKLib/Viterbi.o $(LIBNET) SNetLib/nerror.o SNetLib/nlayer.o SNetLib/nnet.o SNetLib/progobj.o SNetLib/socketobj.o 


###############################################################################
# library specific rules
#


###############################################################################
# cleans the binaries  
clean:
	-rm -f *.o STKLib/*.o SNetLib/*.o SVite SFea2OutP SRover SERest SLRatio SResults SExpand SCombAcc SMMIUpdate SNet
  


###############################################################################
# runs makedepend to automatically define header files dependencies
depend:
	-makedepend -Y *.c *.cc STKLib/*.c STKLib/*.cc SNetLib/*.c SNetLib/*.cc


###############################################################################
# some testing rules
test : $(BASELIBS)
test2 : test2.cc STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/Buffer.o

# DO NOT DELETE

SERest.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SERest.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SERest.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/Viterbi.h
SERest.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/Net.h
SERest.o: STKLib/labels.h STKLib/stkstream.h
SExpand.o: STKLib/Net.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SExpand.o: STKLib/Error.h STKLib/dict.h STKLib/labels.h STKLib/common.h
SFeaCat.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFeaCat.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SFeaCat.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/Viterbi.h
SFeaCat.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/Net.h
SFea2OutP.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFea2OutP.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc
SFea2OutP.o: STKLib/stkstream.h STKLib/stkstream.tcc STKLib/Viterbi.h
SFea2OutP.o: STKLib/Models.h STKLib/labels.h STKLib/dict.h STKLib/Net.h
SLRatio.o: STKLib/Viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SLRatio.o: STKLib/Error.h STKLib/dict.h STKLib/Net.h STKLib/Models.h
SLRatio.o: STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SLRatio.o: STKLib/stkstream.tcc STKLib/fileio.h STKLib/labels.h
SLRatio.o: STKLib/common.h
SResults.o: STKLib/labels.h STKLib/common.h
SRover.o: STKLib/labels.h STKLib/common.h
SVite.o: STKLib/Viterbi.h STKLib/Models.h STKLib/labels.h STKLib/common.h
SVite.o: STKLib/Error.h STKLib/dict.h STKLib/Net.h STKLib/Models.h
SVite.o: STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SVite.o: STKLib/stkstream.tcc STKLib/fileio.h STKLib/labels.h STKLib/common.h
main.o: STKLib/Matrix.h STKLib/Error.h STKLib/Matrix.tcc STKLib/common.h
mainsock.o: SNetLib/socketobj.h STKLib/Error.h
SNet.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SNet.o: STKLib/Models.h STKLib/Matrix.h STKLib/Matrix.tcc STKLib/stkstream.h
SNet.o: STKLib/stkstream.tcc STKLib/Viterbi.h STKLib/Models.h STKLib/labels.h
SNet.o: STKLib/dict.h STKLib/Net.h SNetLib/nnet.h SNetLib/nlayer.h
SNet.o: SNetLib/progobj.h SNetLib/socketobj.h STKLib/Error.h
STKLib/common.o: STKLib/common.h STKLib/Error.h STKLib/filmatch.h
STKLib/dict.o: STKLib/dict.h STKLib/Models.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/stkstream.h STKLib/stkstream.tcc
STKLib/filmatch.o: STKLib/filmatch.h STKLib/common.h STKLib/Error.h
STKLib/getopt.o: STKLib/getopt.h
STKLib/imagesc.o: STKLib/imagesc.h
STKLib/labels.o: STKLib/Viterbi.h STKLib/Models.h STKLib/labels.h
STKLib/labels.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/Net.h
STKLib/Buffer.o: STKLib/Buffer.h
STKLib/Error.o: STKLib/Error.h
STKLib/labelreader.o: STKLib/labelreader.h STKLib/stkstream.h STKLib/common.h
STKLib/labelreader.o: STKLib/Error.h STKLib/stkstream.tcc
STKLib/Matrix.o: STKLib/Matrix.h STKLib/Error.h STKLib/Matrix.tcc
STKLib/Matrix.o: STKLib/common.h
STKLib/Models.o: STKLib/Models.h STKLib/stkstream.h STKLib/common.h
STKLib/Models.o: STKLib/Error.h STKLib/stkstream.tcc
STKLib/Models_IO.o: STKLib/Models.h
STKLib/Net.o: STKLib/Models.h STKLib/labels.h STKLib/common.h STKLib/Error.h
STKLib/Net.o: STKLib/dict.h STKLib/Net.h
STKLib/Net_IO.o: STKLib/labels.h STKLib/common.h STKLib/Error.h STKLib/Net.h
STKLib/Net_IO.o: STKLib/Models.h STKLib/dict.h
STKLib/Viterbi.o: STKLib/Viterbi.h STKLib/Models.h STKLib/labels.h
STKLib/Viterbi.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/Net.h
SNetLib/progobj.o: SNetLib/progobj.h SNetLib/socketobj.h STKLib/Error.h
SNetLib/socketobj.o: SNetLib/socketobj.h STKLib/Error.h

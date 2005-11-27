CFLAGS = -I/usr/X11R6/include -Wall -g -DDOUBLEPRECISION=0
# -DEXPANDNET_ONLY
XFLAGS = -DMOTIF -L/usr/X11R6/lib -lXm -lXt -lX11 -lXpm
#LDFLAGS = -a -g -pg
#ICCFLAGS = -O2 -tpp6 -axiMK -xiMK -parallel
CC  = g++




all  : SVite SRover SERest SResults SExpand SLRatio



SFea2OutP : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/net.o -lm
SVite      : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SLRatio    : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SERest     : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o STKLib/net.o STKLib/dict.o -lm
SFeaCat    : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/fileio.o STKLib/hmms.o -lm
SRover     : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/labels.o -lm
SResults   : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/labels.o -lm
SExpand  : STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/labels.o STKLib/dict.o STKLib/net.o -lm

# STKLib/hmms.o STKLib/viterbi.o

#SCompV     : STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o -lm
#SCombAcc   : STKLib/common.o STKLib/hmms.o -lm
#SMMIUpdate : STKLib/common.o STKLib/hmms.o -lm
#SMMIUpdate_mfe : STKLib/common.o STKLib/hmms.o -lm
#SVite_ctx  : STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o STKLib/net.o -lm
#HMkNet : HTKLib/HTKLib.linux.a -lm
#SCombine :  SCombine.o STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o
#	$(CC) $(CFLAGS) -o SCombine SCombine.o STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.o STKLib/labels.o -lm

#wfst2net   : STKLib/common.o -lm
SERestX :  SERest.o STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.c STKLib/labels.o STKLib/net.o STKLib/dict.o STKLib/imagesc.o
	$(CC) $(CFLAGS) $(XFLAGS) -o SERestX SERest.o STKLib/common.o STKLib/fileio.o STKLib/hmms.o STKLib/viterbi.c STKLib/labels.o STKLib/net.o STKLib/dict.o STKLib/imagesc.o -lm


#STKLib/hmms.o : STKLib/hmms.c STKLib/ModelSet.cc STKLib/Mixture.cc STKLib/Hmm.cc
  
clean:
	-rm -f *.o STKLib/*.o SVite SFea2OutP SRover SERest SLRatio SResults SExpand SCombAcc SMMIUpdate
depend:
	-makedepend -Y *.c HTKLib/*.c STKLib/*.c

test : test.cc STKLib/common.o STKLib/labelreader.o STKLib/filmatch.o STKLib/Buffer.o

test2 : test2.cc STKLib/common.o STKLib/Error.o STKLib/filmatch.o STKLib/Buffer.o


# DO NOT DELETE

SERest.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SERest.o: STKLib/hmms.h STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h
SERest.o: STKLib/dict.h STKLib/net.h STKLib/labels.h
SExpand.o: STKLib/net.h STKLib/hmms.h STKLib/labels.h STKLib/common.h
SExpand.o: STKLib/Error.h STKLib/dict.h STKLib/labels.h STKLib/common.h
SFeaCat.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFeaCat.o: STKLib/hmms.h STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h
SFeaCat.o: STKLib/dict.h STKLib/net.h
SFea2OutP.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h STKLib/common.h
SFea2OutP.o: STKLib/hmms.h STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h
SFea2OutP.o: STKLib/dict.h STKLib/net.h
SLRatio.o: STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h STKLib/common.h
SLRatio.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/hmms.h
SLRatio.o: STKLib/fileio.h STKLib/labels.h STKLib/common.h
SResults.o: STKLib/labels.h STKLib/common.h
SRover.o: STKLib/labels.h STKLib/common.h
SVite.o: STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h STKLib/common.h
SVite.o: STKLib/Error.h STKLib/dict.h STKLib/net.h STKLib/hmms.h
SVite.o: STKLib/fileio.h STKLib/labels.h STKLib/common.h
STKLib/common.o: STKLib/common.h STKLib/Error.h STKLib/filmatch.h
STKLib/dict.o: STKLib/dict.h STKLib/hmms.h STKLib/common.h STKLib/Error.h
STKLib/fileio.o: STKLib/fileio.h STKLib/common.h STKLib/Error.h
STKLib/filmatch.o: STKLib/filmatch.h STKLib/common.h STKLib/Error.h
STKLib/getopt.o: STKLib/getopt.h
STKLib/hmms.o: STKLib/hmms.h STKLib/common.h STKLib/Error.h
STKLib/hmms.o: STKLib/stkstream.h STKLib/stkstream.tcc
STKLib/imagesc.o: STKLib/imagesc.h
STKLib/labels.o: STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h
STKLib/labels.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/net.h
STKLib/net.o: STKLib/hmms.h STKLib/labels.h STKLib/common.h STKLib/Error.h
STKLib/net.o: STKLib/dict.h STKLib/net.h
STKLib/viterbi.o: STKLib/viterbi.h STKLib/hmms.h STKLib/labels.h
STKLib/viterbi.o: STKLib/common.h STKLib/Error.h STKLib/dict.h STKLib/net.h

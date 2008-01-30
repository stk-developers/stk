EXTRA_DIST = doc m4 utils reconf README ToDo.txt 

SUBDIRS = lib src 

dist-hook:
	rm -rf `find $(distdir)/doc -name .svn`

SLRatio SRover SResults SFeaCat SExpand SERest SVite STree SFeaStack:
	cd src; cd STKLib; make; cd ..;  make $@; cd .. 

SNet:
	cd src; cd STKLib; make; cd ../SNetLib; make; cd ..; make $@; cd .. 

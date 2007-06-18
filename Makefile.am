EXTRA_DIST = doc m4 utils reconf README ToDo.txt 

SUBDIRS = lib src 

dist-hook:
	rm -rf `find $(distdir)/doc -name .svn`

STree:
	cd src; make STree; cd .. 

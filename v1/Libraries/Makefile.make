.PHONY : all cleanmodules cleanqrref clean

all:
	(cd Modules;  make )
	(cd libqrref; make )

cleanmodules:
	cd Modules; make clean

cleanqrref:
	cd libqrref; make clean

clean: cleanmodules cleanqrref

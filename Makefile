
all:
	$(MAKE) -C GEMTools
	$(MAKE) -C python

install: all
	$(MAKE) -C python install

install-user: all
	$(MAKE) -C python install-user

clean:
	$(MAKE) -C GEMTools clean
	$(MAKE) -C python clean



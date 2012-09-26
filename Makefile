
all:
	$(MAKE) -C GEMTools
	$(MAKE) -C bindings/python

install: all
	$(MAKE) -C bindings/python install

install-user: all
	$(MAKE) -C bindings/python install-user

clean:
	$(MAKE) -C GEMTools clean
	$(MAKE) -C bindings/python



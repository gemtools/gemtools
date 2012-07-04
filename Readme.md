Gemtools

JIRA : http://barnaserver.com/jira/browse/GT
Wiki : http://barnaserver.com/confluence/x/TAEc
git  : http://barnaserver.com/fisheye/git/gemtools.git


--- make the c library
cd GEMTools
make

-- install the pythin library
cd bindings/python
make
python setup.py install --user

-- gem 2 bed conversion
install the python library
go to bindings/python/examples

cat <your .map file> | ./gem2bed.py > mybed.bed

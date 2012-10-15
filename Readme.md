Gemtools

JIRA : http://barnaserver.com/jira/browse/GT
Wiki : http://barnaserver.com/confluence/x/TAEc
git  : http://barnaserver.com/fisheye/git/gemtools.git

-- Changes --

Release Notes - GEM-Tools - Version 1.2

** Bug
    * [GT-12] - Wrong sam2fastq conversion
    * [GT-14] - The indel parameter (-e) is missing for mapper runs
    * [GT-15] - Merging gem files is messed up due to some wronge reference in the function
    * [GT-8] - Update the merge method to support the new streaming structure
    * [GT-11] - Add the ability to realign SAM files with gem-map-2-map
    * [GT-16] - For pair alignment, the extension should be disabled
    * [GT-18] - Include GEM binaries in the python bundle
    * [GT-7] - Add support to compress files


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




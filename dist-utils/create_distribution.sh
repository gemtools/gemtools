#!/bin/bash

# names and target directories
VERSION=$1
TYPE=$2

NAME=gemtools-$VERSION-$TYPE
DIR=dist/$NAME

# ensure base directlry exists
mkdir -p $DIR
mkdir -p $DIR/lib64
mkdir -p $DIR/bin

if [ ! -e $DIR/lib ]; then
    # create lib link
    $(cd $DIR; ln -s lib64 lib)
fi

# create a virtual environment
# with python 2.6
python dist-utils/virtualenv.py -p python2.6 dist/env-2.6
if [ $? -eq 0 ]; then
    . dist/env-2.6/bin/activate
    # install
    python setup.py clean
    GEM_NO_BUNDLE=1 python setup.py install --prefix=$DIR
    # install argparse for 2.6
    pip install argparse --target=$DIR/lib64/python2.6/site-packages/
    deactivate
fi

# create a virtual environment
# with python 2.7
python dist-utils/virtualenv.py -p python2.7 dist/env-2.7
if [ $? -eq 0 ]; then
    . dist/env-2.7/bin/activate
    # install
    python setup.py clean
    GEM_NO_BUNDLE=1 python setup.py install --prefix=$DIR
    # install argparse for 2.7 just in case
    pip install argparse --target=$DIR/lib64/python2.7/site-packages/
    deactivate
fi

# copy start script
cp -R dist-utils/gemtools.py $DIR/bin/gemtools
chmod +x $DIR/bin/gemtools

# build tools and copy lib/include and bin
make gemtools
cp -R GEMTools/lib/* $DIR/lib64
cp -R GEMTools/include $DIR/include
cp -R GEMTools/bin/* $DIR/bin/

# extract binaries
D26=$DIR/lib64/python2.6/site-packages/gem/gembinaries
D27=$DIR/lib64/python2.7/site-packages/gem/gembinaries
BASE="../../../../.."
if [ -e "$DIR/lib64/python2.6" ]; then mkdir -p $D26; fi
if [ -e "$DIR/lib64/python2.7" ]; then mkdir -p $D27; fi
for f in `tar -C $DIR/bin -xzvf downloads/*-$VERSION-$TYPE*`; do
    echo "Linking $f"
    if [ -e $D26 ]; then rm -f $D26/$f; $(cd $D26; ln -s $BASE/bin/$f .); fi
    if [ -e $D27 ]; then rm -f $D27/$f; $(cd $D27; ln -s $BASE/bin/$f .); fi
done

# create tarball
x=$(cd dist; tar czvf $NAME.tar.gz $NAME)

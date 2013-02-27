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


# copy start script
cp dist/dist/gemtools $DIR/bin/gemtools
chmod +x $DIR/bin/gemtools

# build tools and copy lib/include and bin
make gemtools
cp -R GEMTools/lib/* $DIR/lib64
cp -R GEMTools/include $DIR/include
cp -R GEMTools/bin/* $DIR/bin/

# extract binaries
BASE="../../../../.."
tar -C $DIR/bin -xzvf downloads/*-$VERSION-$TYPE*

# create tarball
x=$(cd dist; tar czvf $NAME.tar.gz $NAME)

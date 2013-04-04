#!/bin/bash
# package the gemtools python library distribution
# 
# We create a dedicatd folder. Install all dependencies to that folder
# and then install gemtools into that folder as well.
# Finally, the binaries are linked to the bin folder and 
# the package is tarred up
#
# Arguments:
#	you have to specify the version number and the type (core2|i3) as
#	arguments
# names and target directories


function download_package {
	cache=$1
	name=$2
	version=$3
	result="${cache}/${name}-${version}.tar.gz"
	if [ ! -e "$result" ]; then
		echo "Downloading ${name}==${version}"
		pip install -d $cache "${name}==${version}"
		if [[ "$?" != "0" ]]; then
			echo "Warning while downloading ${name}-${verison}! If the is matplot lib, ignore this!"
		fi
	fi
}

function prepare_folder {
	DIR=$1
	mkdir -p $DIR
	mkdir -p $DIR/lib64
	mkdir -p $DIR/bin
	if [ ! -e $DIR/lib ]; then
		# create lib link
		$(cd $DIR; ln -s lib64 lib)
	fi
}

function install_dependency {
	cache=$1
	name=$2
	version=$3
	target=$4
	python_version=$5
	s=$6
	result="${cache}/${name}-${version}.tar.gz"
	inst_dir="$target/lib64/python${python_version}/site-packages"
	echo "Installing ${name}-${version}"
	cd $cache && tar xzf $result && cd "${name}-${version}"
   	PYTHONPATH=$PYTHONPATH:${inst_dir} python setup.py install --prefix=$target 
	cd .. && rm -R "${name}-${version}"
	cd $s
}

function install_folder {
	VERSION=$1
	TYPE=$2
	CACHE=$3
	PY=$4

	mkdir -p dist

	# create a virtual environment
	python dist-utils/virtualenv.py -p python${PY} dist/env-${PY}
	if [ $? -eq 0 ]; then
		echo "Preparing folder" 
		NAME=gemtools-$VERSION-$TYPE
		DIR=$(pwd)/dist/$NAME
		prepare_folder $DIR
		. dist/env-$PY/bin/activate
		# install
		python setup.py clean
		GEM_NO_BUNDLE=1 python setup.py install --prefix=$DIR
		INST_DIR=$DIR/lib64/python$PY/site-packages
		CWD=`pwd`
		install_dependency $CACHE  "numpy"  "1.7.0" $DIR $PY $CWD
		install_dependency $CACHE  "argparse"  "1.2.1" $DIR $PY $CWD
		install_dependency $CACHE  "matplotlib"  "1.2.0" $DIR $PY $CWD
		deactivate
	fi
}

VERSION=$1
TYPE=$2
if [[ "$VERSION" == "" ]]; then
	echo "You have to specify the verison"
	exit 1
fi

if [[ "$TYPE" == "" ]]; then
	echo "You have to specify the type as i3 or core2"
	exit 1
fi


# create the cache folder on pwd
CACHE=$(pwd)/downloads
mkdir -p $CACHE

# download the dependencies to the cache
echo "Preparing dependencies"
download_package $CACHE "numpy" "1.7.0"
download_package $CACHE "argparse" "1.2.1"
download_package $CACHE "matplotlib" "1.2.0"

echo "Installing to target"
install_folder $VERSION $TYPE $CACHE "2.6"


## copy start script
cp -R dist-utils/gemtools.py $DIR/bin/gemtools
chmod +x $DIR/bin/gemtools

## build tools and copy lib/include and bin
#make gemtools
cp -R GEMTools/lib/* $DIR/lib64
cp -R GEMTools/include $DIR/include
cp -R GEMTools/bin/* $DIR/bin/

## extract binaries
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

## create tarball
x=$(cd dist; tar czf $NAME.tar.gz $NAME)

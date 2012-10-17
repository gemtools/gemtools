#!/bin/bash
dir=`dirname $0`
home="${dir}/.."
envdir="${home}/gemtools-env"
deactivate 2>/dev/null

rm -Rf $envdir 2>/dev/null

echo "Creating virtual environment"
python $dir/virtualenv.py $envdir
. $envdir/bin/activate

echo "Installing dependencies"
pip install nose

echo "--------------------------------------------"
echo "- Virtual environment created in ${envdir}"
echo "- Load the environment with"
echo "- . ${envdir}/bin/activate"

#!/bin/bash

## Ensure local install
#export CFLAGS="-I$PREFIX/include"
#export LDFLAGS="-L$PREFIX/lib"
#export CPATH=$PREFIX/include

## go into build dir/
cd $SRC_DIR
mkdir -p build/
cd build

sw_vers -productVersion
if [ $? -eq 0 ]
then
    ## see here when writing OSX instructions:
    ## https://github.com/Cibiv/IQ-TREE/wiki/Compilation-Guide#downloading-source-code
    echo "Doing Mac build"
    CC=clang CXX=clang++ /usr/local/bin/cmake ..
else
    ## easy linux install
    echo "Doing linux build"
    cmake -D -DIQTREE_FLAGS=omp ..
fi
make -j4

## pro-tip - I kept trying this with 'conda build .' and afterwards the binary
## was nowhere to be found. So I thought it was failing, even though it said that
## it build successfully. But then I tried conda install -c ipyrad iqtree-pomo
## and it worked. So there you have it. 
cp ./iqtree $PREFIX/bin/iqtree-pomo


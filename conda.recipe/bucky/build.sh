#!/bin/bash

## go to the src dir
cd $SRC_DIR/src

make

cp ./bucky $PREFIX/bin
cp ./mbsum $PREFIX/bin

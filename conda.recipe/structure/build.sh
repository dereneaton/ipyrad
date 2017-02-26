#!/bin/bash

## go to the src dir
cd $SRC_DIR/

## compile
make -f Makefile
rm *.o

## copy binary to bin
cp ./structure $PREFIX/bin

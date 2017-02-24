#!/bin/bash

## go to the src dir
cd $SRC_DIR/src

gcc -o bpp -03 bpp.c tools.c -lm
gcc -o bpp-avx -03 -DUSE_AVX -mavx bpp.c tools.c -lm
gcc -o MCcoal -DSIMULATION bpp.c tools.c -lm 

cp ./bpp $PREFIX/bin
cp ./bpp_avx $PREFIX/bin
cp ./MCcoal $PREFIX/bin

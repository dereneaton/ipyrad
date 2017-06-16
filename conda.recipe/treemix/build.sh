#!/bin/bash

./configure CPPFLAGS=-I${PREFIX}/include LDFLAGS=-L${PREFIX}/lib --prefix=${PREFIX} --with-boost=${PREFIX}
make
make install

cp $SRC_DIR/src/treemix $PREFIX/bin
cp $SRC_DIR/src/threepop $PREFIX/bin
cp $SRC_DIR/src/fourpop $PREFIX/bin
cp $SRC_DIR/src/f4ratio $PREFIX/bin




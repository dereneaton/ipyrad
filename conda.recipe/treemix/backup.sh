#!/bin/bash

./configure CPPFLAGS=-I${PREFIX}/include LDFLAGS=-L${PREFIX}/lib --prefix=${PREFIX} --with-boost=${PREFIX}
make
make install






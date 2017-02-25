#!/bin/bash


$PYTHON setup.py build_ext --inplace
cd $SRC_DIR/vars
$PYTHON setup.py build_ext --inplace

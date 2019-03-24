#!/usr/bin/bash

pushd pop2
pushd build
make
popd
popd
cp pop2/build/libpop.a ../../lib
cp pop2/build/compile/*.mod ../../include

#!/bin/sh
make -j -C ../build/
../build/test/ma_test serialsquare.dmg serialsquare.smb
paraview after/after.pvtu 

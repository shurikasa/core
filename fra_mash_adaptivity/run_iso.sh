#!/bin/sh
make -j -C ../build/
../build/test/ma_test serialsquare.dmg serialsquare.smb 0.1
paraview after/after.pvtu 

#!/bin/sh
make -j -C ../build/
../build/test/aniso_ma_test serialsquare.dmg serialsquare.smb
paraview aniso_after/aniso_after.pvtu 
